NP_Graph = function(W, Z, sigma_eta, rho, sigma_delta, label_name, var_thred = 5) 
{
    correct_W = function(W, sigma_delta) {  ## apply regression calibration under multivariate version
        covariance_matrix_W = as.matrix(stats::cor(W))
        n = nrow(W)
        m = ncol(W)
        Y = matrix(0, nrow = n, ncol = m)
        for (i in 1:n) {
            Y[i, ] = colMeans(W) + t(covariance_matrix_W - sigma_delta) %*% 
                solve(covariance_matrix_W ) %*% 
                (W[i, ] - colMeans(W))
        }
        return(as.data.frame(Y)) }

    correct_Z = function(Z, sigma_eta) {  ## apply regression calibration under multivariate version
        covariance_matrix_Z = as.matrix(stats::cor(Z))
        n = nrow(Z)
        p = ncol(Z)
        X = matrix(0, nrow = n, ncol = p)
        for (i in 1:n) {
            X[i, ] = colMeans(Z) + t(covariance_matrix_Z - sigma_eta) %*% 
                solve(covariance_matrix_Z ) %*% 
                (Z[i, ] - colMeans(Z))
        }
        return(as.data.frame(X))
    }
    Y = as.matrix(correct_W(W, sigma_delta))
    X = as.matrix(correct_Z(Z, sigma_eta))
    models = list()
    importance = list()
    for (i in 1:ncol(Y)) {
        df_Y = Y[, i]
        df = cbind(X, Y = df_Y)
        model = randomForest::randomForest(Y ~ ., data = df, 
            importance = TRUE)
        models[[i]] = model
        importance[[i]] = caret::varImp(model)
    }
    importance_matrix = do.call(cbind, importance)
    rf_model_refit = list()
    for (i in 1:ncol(Y)) {
        df_Y = Y[, i]
        df = cbind(X, Y = df_Y)
        imp_sel = importance[[i]][, "Overall"] > var_thred
        important_vars = rownames(importance[[i]][imp_sel, , 
            drop = FALSE])
        if (length(important_vars) > 0) {
            rf_model_refit[[i]] = randomForest::randomForest(Y ~ 
                ., data = df[, c("Y", important_vars)])
        }
        else {
            max_var = rownames(importance[[i]])[which.max(importance[[i]][, 
                "Overall"])]
            rf_model_refit[[i]] = randomForest::randomForest(Y ~ 
                ., data = df[, c("Y", max_var)])
        }
    }
    r_list = list()
    for (i in 1:ncol(Y)) {
        r_i = Y[, i] - stats::predict(rf_model_refit[[i]], X)
        r_list[[i]] = r_i
    }
    r_matrix = do.call(cbind, r_list)
    PSE = mean(sapply(r_list, function(r_matrix) sum(r_matrix^2)))
    y_matrix = as.matrix(Y)
    n_resp = ncol(y_matrix)
    correlation = matrix(1, nrow = n_resp, ncol = n_resp)
    for (i in 1:n_resp) {
        for (j in 1:n_resp) {
            if (i != j) {
                correlation[i, j] = metrica::dcorr(obs = r_matrix[, 
                  i], pred = r_matrix[, j])
            }
        }
    }
    for (i in 1:(n_resp - 1)) {
        for (j in (i + 1):n_resp) {
            max_value = max(correlation[i, j], correlation[j, 
                i])
            correlation[i, j] = max_value
            correlation[j, i] = max_value
        }
    }
    glasso_result = glasso::glasso(correlation, rho = rho)
    precision_matrix = glasso_result$wi
    net = precision_matrix
    net = network::network(net, directed = FALSE)
    network::network.vertex.names(net) = paste0("Y", network::network.vertex.names(net))
    graph = GGally::ggnet2(net, size = 10, node.color = "lightgray", 
        label = label_name, label.size = 3, mode = "circle")
    return(list(W_hat = Y, Z_hat = X, PSE = PSE, importance_matrix = importance_matrix, 
        precision_matrix = precision_matrix, graph = graph))
}
