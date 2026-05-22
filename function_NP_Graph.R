NP_Graph = function(W, Z, sigma_eta, rho, sigma_delta = 0.5, r = 0.8, lambda = 1, 
    pi = 0.8, label_name, var_thred = 5) 
{
    correct_W = function(W, sigma_delta, r, lambda, pi) {
        m = ncol(W)
        n = nrow(W)
        model = c()
        for (j in 1:m) {
            if (sum(W[, j]%%1 == 0) == n & sum(W[, j] != 0 & 
                W[, j] != 1) > 0) {
                model = c(model, "counts")
            }
            else if (sum(W[, j]%%1 == 0) == n & sum(W[, j] != 
                0 & W[, j] != 1) == 0) {
                model = c(model, "binary")
            }
            else {
                model = c(model, "continuous")
            }
        }
        Y = data.frame(matrix(ncol = m, nrow = n))
        colnames(Y) = colnames(W)
        for (j in 1:m) {
            if (model[j] == "continuous") {
                Y[, j] = mean(W[, j]) + (stats::var(W[, j]) - 
                  sigma_delta[j,j])/stats::var(W[, j]) * (W[, j] - 
                  mean(W[, j]))
            }
            else if (model[j] == "binary") {
                S = stats::rbinom(n, 1, r)
                Y[, j] = (W[, j] + S - 1)/(2 * S - 1)
            }
            else {
                Y[, j] = ((mean(W[, j]) - lambda)/(1 - pi)) + 
                  ((mean(W[, j]) - lambda)/(1 - pi)) * (mean(W[, 
                    j]) - lambda)/(lambda + (3 * pi + 1)/(1 - 
                    pi) * (mean(W[, j]) - lambda)) * (W[, j] - 
                    mean(W[, j]))
            }
        }
        return(Y)
    }
    correct_Z = function(Z, sigma_eta) {
        covariance_matrix_Z = as.matrix(stats::cor(Z))
        n = nrow(Z)
        p = ncol(Z)
        X = matrix(0, nrow = n, ncol = p)
        for (i in 1:n) {
            X[i, ] = colMeans(Z) + t(covariance_matrix_Z - sigma_eta) %*% 
                solve(covariance_matrix_Z + diag(0.2, dim(Z)[2])) %*% 
                (Z[i, ] - colMeans(Z))
        }
        return(as.data.frame(X))
    }
    Y = as.matrix(correct_W(W, sigma_delta, r, lambda, pi))
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
                correlation[i, j] = metrica::dcorr(obs = y_matrix[, 
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
