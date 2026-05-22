#library(PAGE)

## function setup

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


##############################################################################
## Analysis Start
library(lattice)
library(randomForest)
library(ggplot2)

Y = read.table("C://final_microRNA.csv",sep=",",head=TRUE)
Y = Y[,-1]
mRNA = colnames(Y)

X = read.table("C://final_GeneExpression.csv",sep=",",head=TRUE)
X = X[,-1]
gene = colnames(X)

n = dim(Y)[1]
m = dim(Y)[2]
p = dim(X)[2]

Y = matrix(unlist(Y),n,m)
X = matrix(unlist(X),n,p)

### Sensitivity analyses

R = 0.5; R1 = 0.8
Sigma_Y = cov(Y)
Sigma_X = cov(X)

DA1 = NP_Graph(Y, X,sigma_eta=(1-R)*Sigma_X,rho=.2,sigma_delta = (1-R)*Sigma_Y, label=mRNA)
DA2 = NP_Graph(Y, X,sigma_eta=(1-R1)*Sigma_X,rho=.2,sigma_delta = (1-R1)*Sigma_Y, label=mRNA)
DA3 = NP_Graph(Y, X,sigma_eta=diag(0.2,p),rho=.2,sigma_delta = diag(0.2,m), label=mRNA)
DA4 = NP_Graph(Y, X,sigma_eta=diag(0.5,p),rho=.2,sigma_delta = diag(0.5,m), label=mRNA)
DA5 = NP_Graph(Y, X,sigma_eta=diag(0,p),rho=.2,sigma_delta = diag(0,m), label=mRNA)

M1 = DA1$precision_matrix
M2 = DA2$precision_matrix
M3 = DA3$precision_matrix
M4 = DA4$precision_matrix
M5 = DA5$precision_matrix

M1[which(abs(M1)<0.1)]=0
M2[which(abs(M2)<0.1)]=0
M3[which(abs(M3)<0.1)]=0
M4[which(abs(M4)<0.1)]=0
M5[which(abs(M5)<0.1)]=0


##############################
## detection of network structures

  net = M1
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("Y",network::network.vertex.names(net))
  GGally::ggnet2(net,size=10,node.color = "lightgray",label=mRNA,label.size = 2.7,mode = "circle")+ 
  labs(title = "Scenario II-0.5")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


  net = M2
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("Y",network::network.vertex.names(net))
  GGally::ggnet2(net,size=10,node.color = "lightgray",label=mRNA,label.size = 2.7,mode = "circle")+ 
  labs(title = "Scenario II-0.8")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  net = M3
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("Y",network::network.vertex.names(net))
  GGally::ggnet2(net,size=10,node.color = "lightgray",label=mRNA,label.size = 2.7,mode = "circle")+ 
  labs(title = "Scenario I-0.2")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


  net = M4
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("Y",network::network.vertex.names(net))
  GGally::ggnet2(net,size=10,node.color = "lightgray",label=mRNA,label.size = 2.7,mode = "circle")+ 
  labs(title = "Scenario I-0.5")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  net = M5
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("Y",network::network.vertex.names(net))
  GGally::ggnet2(net,size=10,node.color = "lightgray",label=mRNA,label.size = 2.7,mode = "circle")+ 
  labs(title = "naive")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

##############################
# identification of the most important variable

B1 = c(); B2 = c(); B3 = c(); B4 = c(); B5 = c()
id1 = c(); id2 = c(); id3 = c(); id4 = c(); id5 = c()
est1 = DA1$importance_matrix
est2 = DA2$importance_matrix
est3 = DA3$importance_matrix
est4 = DA4$importance_matrix
est5 = DA5$importance_matrix

for(i in 1:20){

id1 = c(id1,which(est1[,i] == max(est1[,i])))
id2 = c(id2,which(est2[,i] == max(est2[,i])))
id3 = c(id3,which(est3[,i] == max(est3[,i])))
id4 = c(id4,which(est4[,i] == max(est4[,i])))
id5 = c(id5,which(est5[,i] == max(est5[,i])))

B1 = c(B1,gene[which(est1[,i] == max(est1[,i]))])
B2 = c(B2,gene[which(est2[,i] == max(est2[,i]))])
B3 = c(B3,gene[which(est3[,i] == max(est3[,i]))])
B4 = c(B4,gene[which(est4[,i] == max(est4[,i]))])
B5 = c(B5,gene[which(est5[,i] == max(est5[,i]))])

}

#B = rbind(B1,B2)
#colnames(B) = mRNA
#barchart(as.vector(B))
barchart(as.vector(B1),main="Scenario II-0.5")
barchart(as.vector(B2),main="Scenario II-0.8")
barchart(as.vector(B3),main="Scenario I-0.2")
barchart(as.vector(B4),main="Scenario I-0.5")
barchart(as.vector(B5),main="naive")


##############################################
## compute predicted square errors


PSE1 = 0; PSE2 = 0; PSE3 = 0; PSE4 = 0; PSE5 = 0

IM1 = matrix(unlist(est1),500,20); IM1[which(abs(IM1)<4)] = 0; length(which(IM1!=0))
IM2 = matrix(unlist(est2),500,20); IM2[which(abs(IM2)<4)] = 0; length(which(IM2!=0))
IM3 = matrix(unlist(est3),500,20); IM3[which(abs(IM3)<4)] = 0; length(which(IM3!=0))
IM4 = matrix(unlist(est4),500,20); IM4[which(abs(IM4)<4)] = 0; length(which(IM4!=0))
IM5 = matrix(unlist(est5),500,20); IM5[which(abs(IM5)<4)] = 0; length(which(IM5!=0))


sigma_etaR = (1-R)*Sigma_X; sigma_deltaR = (1-R)*Sigma_Y
Z = NULL
for(i in 1:20) {
W =  mean(Y[, i]) + (var(Y[, i]) - sigma_deltaR[i,i]) / var(Y[, i]) * (Y[, i] - mean(Y[, i]))
col = IM1[,i]
for(j in 1:500) {
if(col[j]!=0)
Z = cbind(Z, mean(X[, j]) + (var(X[, j]) - sigma_etaR[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j])))
}

DATA = data.frame(W, Z)
model = randomForest(W ~ ., data = DATA)
pred = predict(model)

PSE1 = PSE1 + sum(abs(pred - W))

}

#########
sigma_etaR1 = (1-R1)*Sigma_X; sigma_deltaR1 = (1-R1)*Sigma_Y
Z = NULL
for(i in 1:20) {
W =  mean(Y[, i]) + (var(Y[, i]) - sigma_deltaR1[i,i]) / var(Y[, i]) * (Y[, i] - mean(Y[, i]))
col = IM2[,i]
for(j in 1:500) {
if(col[j]!=0)
Z = cbind(Z, mean(X[, j]) + (var(X[, j]) - sigma_etaR1[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j])))
}

DATA = data.frame(W, Z)
model = randomForest(W ~ ., data = DATA)
pred = predict(model)

PSE2 = PSE2 + sum(abs(pred - W))

}

#########
Z = NULL
for(i in 2:20) {
W =  mean(Y[, i]) + (var(Y[, i]) - 0.2) / var(Y[, i]) * (Y[, i] - mean(Y[, i]))
col = IM3[,i]
for(j in 1:500) {
if(col[j]!=0)
Z = cbind(Z, mean(X[, j]) + (var(X[, j]) - 0.2) / var(X[, j]) * (X[, j] - mean(X[, j])))
}

DATA = data.frame(W, Z)
model = randomForest(W ~ ., data = DATA)
pred = predict(model)

PSE3 = PSE3 + sum(abs(pred - W))

}

#########
Z = NULL
for(i in 1:20) {
W =  mean(Y[, i]) + (var(Y[, i]) - 0.5) / var(Y[, i]) * (Y[, i] - mean(Y[, i]))
col = IM4[,i]
for(j in 1:500) {
if(col[j]!=0)
Z = cbind(Z, mean(X[, j]) + (var(X[, j]) - 0.5) / var(X[, j]) * (X[, j] - mean(X[, j])))
}

DATA = data.frame(W, Z)
model = randomForest(W ~ ., data = DATA)
pred = predict(model)

PSE4 = PSE4 + sum(abs(pred - W))

}

#########
Z = NULL
for(i in 1:20) {
W =  mean(Y[, i]) + (var(Y[, i]) - 0) / var(Y[, i]) * (Y[, i] - mean(Y[, i]))
col = IM5[,i]
for(j in 1:500) {
if(col[j]!=0)
Z = cbind(Z, mean(X[, j]) + (var(X[, j]) - 0) / var(X[, j]) * (X[, j] - mean(X[, j])))
}

DATA = data.frame(W, Z)
model = randomForest(W ~ ., data = DATA)
pred = predict(model)

PSE5 = PSE5 + sum(abs(pred - W))

}

round(PSE1 / (dim(Y)[1]*dim(Y)[2]),3)
round(PSE2 / (dim(Y)[1]*dim(Y)[2]),3)
round(PSE3 / (dim(Y)[1]*dim(Y)[2]),3)
round(PSE4 / (dim(Y)[1]*dim(Y)[2]),3)
round(PSE5 / (dim(Y)[1]*dim(Y)[2]),3)

#############################################################
##     estimate the nonlinear curves-microRNA vs genes     ##
#############################################################

##Scenario I-0.2
smooth_data_I_02 = NULL
j = which(gene == "GABBR2")
for(i in 1:length(which(B3 == "GABBR2"))) {
k = which(B3 == "GABBR2")[i]
Z =  mean(X[, j]) + (var(X[, j]) - 0.2) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - 0.2) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_I_02[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_I_02[[1]]$x, smooth_data_I_02[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B3 == "GABBR2")[1]]), lwd = 2, ylim = c(-0.8,0.6), xlim=c(-2,4))

plot(smooth_data_I_02[[2]]$x, smooth_data_I_02[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B3 == "GABBR2")[2]]), lwd = 2, ylim = c(-1.3,0), xlim=c(-2,4))


plot(smooth_data_I_02[[3]]$x, smooth_data_I_02[[3]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B3 == "GABBR2")[3]]), lwd = 2, ylim = c(-0.6,0.3), xlim=c(-2,3))


###

smooth_data_I_02_2 = NULL
j = which(gene == "TF")
for(i in 1:length(which(B3 == "TF"))) {
k = which(B3 == "TF")[i]
Z =  mean(X[, j]) + (var(X[, j]) - 0.2) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - 0.2) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_I_02_2[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_I_02_2[[1]]$x, smooth_data_I_02_2[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "TF",
     ylab = paste(mRNA[which(B3 == "TF")[1]]), lwd = 2)

plot(smooth_data_I_02_2[[2]]$x, smooth_data_I_02_2[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "TF",
     ylab = paste(mRNA[which(B3 == "TF")[2]]), lwd = 2)



##########################

##Scenario I-0.5
smooth_data_I_05 = NULL
j = which(gene == "GABBR2")
for(i in 1:length(which(B4 == "GABBR2"))) {
k = which(B4 == "GABBR2")[i]
Z =  mean(X[, j]) + (var(X[, j]) - 0.5) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - 0.5) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_I_05[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_I_05[[1]]$x, smooth_data_I_05[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B4 == "GABBR2")[1]]), lwd = 2, ylim = c(-0.8,0.6), xlim=c(-2,4))

plot(smooth_data_I_05[[2]]$x, smooth_data_I_05[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B4 == "GABBR2")[2]]), lwd = 2, ylim = c(-1.3,0), xlim=c(-2,4))


plot(smooth_data_I_05[[3]]$x, smooth_data_I_05[[3]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B4 == "GABBR2")[3]]), lwd = 2, ylim = c(-0.6,0.3), xlim=c(-2,3))


##########################

##Scenario II-0.5
smooth_data_II_05 = NULL
j = which(gene == "HOXA4")
for(i in 1:length(which(B1 == "HOXA4"))) {
k = which(B1 == "HOXA4")[i]
Z =  mean(X[, j]) + (var(X[, j]) - sigma_etaR[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - sigma_deltaR[k,k]) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_05[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_05[[1]]$x, smooth_data_II_05[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "HOXA4",
     ylab = paste(mRNA[which(B1 == "HOXA4")[1]]), lwd = 2)

plot(smooth_data_II_05[[2]]$x, smooth_data_II_05[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "HOXA4",
     ylab = paste(mRNA[which(B1 == "HOXA4")[2]]), lwd = 2)



#####

smooth_data_II_05_2 = NULL
j = which(gene == "C20orf103")
for(i in 1:length(which(B1 == "C20orf103"))) {
k = which(B1 == "C20orf103")[i]
Z =  mean(X[, j]) + (var(X[, j]) - sigma_etaR[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - sigma_deltaR[k,k]) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_05_2[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_05_2[[1]]$x, smooth_data_II_05_2[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "C20orf103",
     ylab = paste(mRNA[which(B1 == "C20orf103")[1]]), lwd = 2)

plot(smooth_data_II_05_2[[2]]$x, smooth_data_II_05_2[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "C20orf103",
     ylab = paste(mRNA[which(B1 == "C20orf103")[2]]), lwd = 2)



##########################

##Scenario II-0.8
smooth_data_II_08 = NULL
j = which(gene == "VSNL1")
for(i in 1:length(which(B2 == "VSNL1"))) {
k = which(B2 == "VSNL1")[i]
Z =  mean(X[, j]) + (var(X[, j]) - sigma_etaR1[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - sigma_deltaR1[k,k]) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_08[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_08[[1]]$x, smooth_data_II_08[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "VSNL1",
     ylab = paste(mRNA[which(B2 == "VSNL1")[1]]), lwd = 2)

plot(smooth_data_II_08[[2]]$x, smooth_data_II_08[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "VSNL1",
     ylab = paste(mRNA[which(B2 == "VSNL1")[2]]), lwd = 2)



#####

smooth_data_II_08_02 = NULL
j = which(gene == "RND3")
for(i in 1:length(which(B2 == "RND3"))) {
k = which(B2 == "RND3")[i]
Z =  mean(X[, j]) + (var(X[, j]) - sigma_etaR1[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - sigma_deltaR1[k,k]) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_08_02[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_08_02[[1]]$x, smooth_data_II_08_02[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "RND3",
     ylab = paste(mRNA[which(B2 == "RND3")[1]]), lwd = 2)

plot(smooth_data_II_08_02[[2]]$x, smooth_data_II_08_02[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "RND3",
     ylab = paste(mRNA[which(B2 == "RND3")[2]]), lwd = 2)


#####

smooth_data_II_08_03 = NULL
j = which(gene == "LPL")
for(i in 1:length(which(B2 == "LPL"))) {
k = which(B2 == "LPL")[i]
Z =  mean(X[, j]) + (var(X[, j]) - sigma_etaR1[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - sigma_deltaR1[k,k]) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_08_03[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_08_03[[1]]$x, smooth_data_II_08_03[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "LPL",
     ylab = paste(mRNA[which(B2 == "LPL")[1]]), lwd = 2)

plot(smooth_data_II_08_03[[2]]$x, smooth_data_II_08_03[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "LPL",
     ylab = paste(mRNA[which(B2 == "LPL")[2]]), lwd = 2)


#####

smooth_data_II_08_04 = NULL
j = which(gene == "HLA.DPB1")
for(i in 1:length(which(B2 == "LPL"))) {
k = which(B2 == "HLA.DPB1")[i]
Z =  mean(X[, j]) + (var(X[, j]) - sigma_etaR1[j,j]) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - sigma_deltaR1[k,k]) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_08_04[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_08_04[[1]]$x, smooth_data_II_08_04[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "HLA.DPB1",
     ylab = paste(mRNA[which(B2 == "HLA.DPB1")[1]]), lwd = 2)

plot(smooth_data_II_08_04[[2]]$x, smooth_data_II_08_04[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "HLA.DPB1",
     ylab = paste(mRNA[which(B2 == "HLA.DPB1")[2]]), lwd = 2)






##Naive
smooth_data_naive = NULL
j = which(gene == "GABBR2")
for(i in 1:length(which(B5 == "GABBR2"))) {
k = which(B5 == "GABBR2")[i]
Z =  mean(X[, j]) + (var(X[, j]) - 0) / var(X[, j]) * (X[, j] - mean(X[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - 0) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_naive[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_naive[[1]]$x, smooth_data_naive[[1]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B5 == "GABBR2")[1]]), lwd = 2, ylim = c(-0.8,0.6), xlim=c(-2,4))

plot(smooth_data_naive[[2]]$x, smooth_data_naive[[2]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B5 == "GABBR2")[2]]), lwd = 2, ylim = c(-1.3,0), xlim=c(-2,4))


plot(smooth_data_I_05[[3]]$x, smooth_data_I_05[[3]]$y,
     type = "l",      # "l" 代表畫線
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B5 == "GABBR2")[3]]), lwd = 2, ylim = c(-0.6,0.3), xlim=c(-2,3))





################################################################
##     estimate the nonlinear curves-microRNA vs microRNA     ##
################################################################


microRNA = function(j,k,Sigma_delta) {

Z =  mean(Y[, j]) + (var(Y[, j]) - Sigma_delta[j,j]) / var(Y[, j]) * (Y[, j] - mean(Y[, j]))
W =  mean(Y[, k]) + (var(Y[, k]) - Sigma_delta[k,k]) / var(Y[, k]) * (Y[, k] - mean(Y[, k]))
DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_microRNA = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

return(smooth_data_microRNA)

}

##################################### Uniqueness by Scenario II

j = which(mRNA == "hsa.mir.204")
k = which(mRNA == "hsa.mir.222")
smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR1)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,3), ylim=c(-0.8,1))

smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,3), ylim=c(-0.8,1))


##################################### Common triangle

j = which(mRNA == "hsa.mir.136")
k = which(mRNA == "hsa.mir.377")
smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR1)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-1.8,1.8))

smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-1.8,1.8))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.2,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-1.8,1.8))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.5,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-1.8,1.8))

smooth_data = microRNA(j,k,Sigma_delta = diag(0,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-1.8,1.8))

######

j = which(mRNA == "hsa.mir.377")
k = which(mRNA == "hsa.mir.376a")
smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR1)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.5,2.5))

smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.5,2.5))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.2,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.5,2.5))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.5,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.5,2.5))

smooth_data = microRNA(j,k,Sigma_delta = diag(0,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.5,2.5))


######

j = which(mRNA == "hsa.mir.136")
k = which(mRNA == "hsa.mir.376a")
smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR1)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.1,1.8))

smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.1,1.8))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.2,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.1,1.8))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.5,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.1,1.8))

smooth_data = microRNA(j,k,Sigma_delta = diag(0,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,2.8), ylim=c(-2.1,1.8))




##################################### Uniqueness by Scenario I

j = which(mRNA == "hsa.mir.630")
k = which(mRNA == "hsa.mir.7")
smooth_data = microRNA(j,k,Sigma_delta = diag(0.2,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,4.5), ylim=c(-0.8,1))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.5,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,4.5), ylim=c(-0.8,1))

smooth_data = microRNA(j,k,Sigma_delta = diag(0,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,4.5), ylim=c(-0.8,1))

######
j = which(mRNA == "hsa.mir.10b")
k = which(mRNA == "hsa.mir.222")

smooth_data = microRNA(j,k,Sigma_delta = diag(0.2,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,3), ylim=c(-1.8,1))

smooth_data = microRNA(j,k,Sigma_delta = diag(0.5,m))
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,3), ylim=c(-1.8,1))

smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR1)
plot(smooth_data$x, smooth_data$y,
     type = "l",      # "l" 代表畫線
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,3), ylim=c(-1.8,1))








ord_1 = order(Z)
smooth_data_1 = lowess(Z[ord_1], pred[ord_1], f = 0.5) 
# 畫線圖
plot(smooth_data_1$x, smooth_data_1$y,
     type = "l",      # "l" 代表畫線
     xlab = "hsa.mir.376a",  #k=14
     ylab = "hsa.mir.377",  #j=15
     lwd = 2#,
#     ylim = c(-1,3.2)
)














