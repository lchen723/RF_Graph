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
col = IM4[,i]
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


correct_data = function(W, Z, sigma_eta, sigma_delta) 
{
    correct_W = function(W, sigma_delta) {
        covariance_matrix_W = as.matrix(stats::cor(W))
        n = nrow(W)
        m = ncol(W)
        Y = matrix(0, nrow = n, ncol = m)
        for (i in 1:n) {
            Y[i, ] = colMeans(W) + t(covariance_matrix_W - sigma_delta) %*% 
                solve(covariance_matrix_W + diag(0.02, dim(W)[2])) %*% 
                (W[i, ] - colMeans(W))
        }
        return(as.data.frame(Y)) }

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
Y = correct_W(W, sigma_delta)
X = correct_W(Z, sigma_eta)
Data = list(Y,X); names(Data) = c("Y","X")
return(Data)

}

###################################

##Scenario I-0.2
I_02 = correct_data(Y,X,sigma_eta=diag(0.2,p),sigma_delta = diag(0.2,m))

smooth_data_I_02 = NULL
j = which(gene == "GABBR2")
for(i in 1:length(which(B3 == "GABBR2"))) {
k = which(B3 == "GABBR2")[i]
Z =  I_02$X[,j]
W =  I_02$Y[,k]


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_I_02[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_I_02[[1]]$x, smooth_data_I_02[[1]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B3 == "GABBR2")[1]]), lwd = 2, ylim = c(-1,1.6), xlim=c(-13,13))

plot(smooth_data_I_02[[2]]$x, smooth_data_I_02[[2]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B3 == "GABBR2")[2]]), lwd = 2, ylim = c(-1.2,2), xlim=c(-13,13))


###

smooth_data_I_02_2 = NULL
j = which(gene == "RPL39L")
for(i in 1:length(which(B3 == "RPL39L"))) {
k = which(B3 == "RPL39L")[i]
Z =  I_02$X[,j]
W =  I_02$Y[,k]

DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_I_02_2[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_I_02_2[[1]]$x, smooth_data_I_02_2[[1]]$y,
     type = "l",      
     xlab = "RPL39L",
     ylab = paste(mRNA[which(B3 == "RPL39L")[1]]), lwd = 2)

plot(smooth_data_I_02_2[[2]]$x, smooth_data_I_02_2[[2]]$y,
     type = "l",      
     xlab = "RPL39L",
     ylab = paste(mRNA[which(B3 == "RPL39L")[2]]), lwd = 2)



##########################

##Scenario I-0.5
I_05 = correct_data(Y,X,sigma_eta=diag(0.5,p),sigma_delta = diag(0.5,m))

smooth_data_I_05 = NULL
j = which(gene == "GABBR2")
for(i in 1:length(c(which(B4 == "GABBR2"),15))) {
k = c(which(B4 == "GABBR2"),15)[i]
Z =  I_05$X[,j]
W =  I_05$Y[,k]


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_I_05[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_I_05[[1]]$x, smooth_data_I_05[[1]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B4 == "GABBR2")[1]]), lwd = 2, ylim = c(-1,1.6), xlim=c(-13,13))

plot(smooth_data_I_05[[2]]$x, smooth_data_I_05[[2]]$y,
     type = "l",     
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B4 == "GABBR2")[2]]), lwd = 2)


plot(smooth_data_I_05[[3]]$x, smooth_data_I_05[[3]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B4 == "GABBR2")[3]]), lwd = 2)

plot(smooth_data_I_05[[4]]$x, smooth_data_I_05[[4]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B4 == "GABBR2")[4]]), lwd = 2)

plot(smooth_data_I_05[[5]]$x, smooth_data_I_05[[5]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[15]), lwd = 2, ylim = c(-1.2,2), xlim=c(-13,13))


##########################

##Scenario II-0.5
II_05 = correct_data(Y,X,sigma_eta=(1-R)*Sigma_X, sigma_delta = (1-R)*Sigma_Y)

smooth_data_II_05 = NULL
j = which(gene == "VSNL1")
for(i in 1:length(which(B1 == "VSNL1"))) {
k = which(B1 == "VSNL1")[i]
Z =  II_05$X[,j]
W =  II_05$Y[,k]


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_05[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_05[[1]]$x, smooth_data_II_05[[1]]$y,
     type = "l",      
     xlab = "VSNL1",
     ylab = paste(mRNA[which(B1 == "VSNL1")[1]]), lwd = 2)

plot(smooth_data_II_05[[2]]$x, smooth_data_II_05[[2]]$y,
     type = "l",      
     xlab = "VSNL1",
     ylab = paste(mRNA[which(B1 == "VSNL1")[2]]), lwd = 2)



#####

smooth_data_II_05_2 = NULL
j = which(gene == "TF")
for(i in 1:length(which(B1 == "TF"))) {
k = which(B1 == "TF")[i]
Z =  II_05$X[,j]
W =  II_05$Y[,k]

DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_05_2[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_05_2[[1]]$x, smooth_data_II_05_2[[1]]$y,
     type = "l",      
     xlab = "TF",
     ylab = paste(mRNA[which(B1 == "TF")[1]]), lwd = 2)

plot(smooth_data_II_05_2[[2]]$x, smooth_data_II_05_2[[2]]$y,
     type = "l",     
     xlab = "TF",
     ylab = paste(mRNA[which(B1 == "TF")[2]]), lwd = 2)


plot(smooth_data_II_05_2[[3]]$x, smooth_data_II_05_2[[3]]$y,
     type = "l",      
     xlab = "TF",
     ylab = paste(mRNA[which(B1 == "TF")[3]]), lwd = 2)

####

smooth_data_II_05_3 = NULL
j = which(gene == "MAL")
for(i in 1:length(which(B1 == "MAL"))) {
k = which(B1 == "MAL")[i]
Z =  II_05$X[,j]
W =  II_05$Y[,k]

DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_05_3[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_05_3[[1]]$x, smooth_data_II_05_3[[1]]$y,
     type = "l",      
     xlab = "MAL",
     ylab = paste(mRNA[which(B1 == "MAL")[1]]), lwd = 2)

plot(smooth_data_II_05_3[[2]]$x, smooth_data_II_05_3[[2]]$y,
     type = "l",      
     xlab = "MAL",
     ylab = paste(mRNA[which(B1 == "MAL")[2]]), lwd = 2)




##########################

##Scenario II-0.8
II_08 = correct_data(Y,X,sigma_eta=(1-R1)*Sigma_X, sigma_delta = (1-R1)*Sigma_Y)

smooth_data_II_08 = NULL
j = which(gene == "LEFTY2")
for(i in 1:length(which(B2 == "LEFTY2"))) {
k = which(B2 == "LEFTY2")[i]
Z =  II_08$X[,j]
W =  II_08$Y[,k]

DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_08[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_08[[1]]$x, smooth_data_II_08[[1]]$y,
     type = "l",      
     xlab = "LEFTY2",
     ylab = paste(mRNA[which(B2 == "LEFTY2")[1]]), lwd = 2)

plot(smooth_data_II_08[[2]]$x, smooth_data_II_08[[2]]$y,
     type = "l",      
     xlab = "LEFTY2",
     ylab = paste(mRNA[which(B2 == "LEFTY2")[2]]), lwd = 2)



#####

smooth_data_II_08_02 = NULL
j = which(gene == "BBOX1")
for(i in 1:length(which(B2 == "BBOX1"))) {
k = which(B2 == "BBOX1")[i]
Z =  II_08$X[,j]
W =  II_08$Y[,k]

DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_II_08_02[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_II_08_02[[1]]$x, smooth_data_II_08_02[[1]]$y,
     type = "l",      
     xlab = "BBOX1",
     ylab = paste(mRNA[which(B2 == "BBOX1")[1]]), lwd = 2)

plot(smooth_data_II_08_02[[2]]$x, smooth_data_II_08_02[[2]]$y,
     type = "l",      
     xlab = "BBOX1",
     ylab = paste(mRNA[which(B2 == "BBOX1")[2]]), lwd = 2)




##Naive
Naive = correct_data(Y,X,sigma_eta=diag(0,p),sigma_delta = diag(0,m))

smooth_data_naive = NULL
j = which(gene == "GABBR2")
for(i in 1:length(which(B5 == "GABBR2"))) {
k = which(B5 == "GABBR2")[i]
Z =  Naive$X[,j]
W =  Naive$Y[,k]

DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_naive[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_naive[[1]]$x, smooth_data_naive[[1]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B5 == "GABBR2")[1]]), lwd = 2, ylim = c(-1,1.6), xlim=c(-13,13))

plot(smooth_data_naive[[2]]$x, smooth_data_naive[[2]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B5 == "GABBR2")[2]]), lwd = 2)


plot(smooth_data_I_05[[3]]$x, smooth_data_I_05[[3]]$y,
     type = "l",      
     xlab = "GABBR2",
     ylab = paste(mRNA[which(B5 == "GABBR2")[3]]), lwd = 2, ylim = c(-1.2,2), xlim=c(-13,13))



smooth_data_naive = NULL
j = which(gene == "CHGA")
for(i in 1:length(which(B5 == "CHGA"))) {
k = which(B5 == "CHGA")[i]
Z =  Naive$X[,j]
W =  Naive$Y[,k]


DATA = data.frame(W, Z)
model = randomForest(W ~ Z, data = DATA)
pred = predict(model)

ord_1 = order(Z)
smooth_data_naive[[i]] = lowess(Z[ord_1], pred[ord_1], f = 0.5) 

}

plot(smooth_data_naive[[1]]$x, smooth_data_naive[[1]]$y,
     type = "l",     
     xlab = "CHGA",
     ylab = paste(mRNA[which(B5 == "CHGA")[1]]), lwd = 2)

plot(smooth_data_naive[[2]]$x, smooth_data_naive[[2]]$y,
     type = "l",     
     xlab = "CHGA",
     ylab = paste(mRNA[which(B5 == "CHGA")[2]]), lwd = 2)











################################################################
##     estimate the nonlinear curves-microRNA vs microRNA     ##
################################################################


microRNA = function(W,Z) {

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
     type = "l",     
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,3), ylim=c(-0.8,1))

smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR)
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-2,3), ylim=c(-0.8,1))


##################################### Common triangle

j = which(mRNA == "hsa.mir.136")
k = which(mRNA == "hsa.mir.377")
smooth_data = microRNA(II_08$Y[,k], II_08$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,5), ylim=c(-4,3))

smooth_data = microRNA(II_05$Y[,k], II_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,5), ylim=c(-4,3))

smooth_data = microRNA(I_05$Y[,k], I_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,5), ylim=c(-4,3))

smooth_data = microRNA(I_02$Y[,k], I_02$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,5), ylim=c(-4,3))

smooth_data = microRNA(Naive$Y[,k], Naive$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,5), ylim=c(-4,3))

######

j = which(mRNA == "hsa.mir.377")
k = which(mRNA == "hsa.mir.376a")
smooth_data = microRNA(II_08$Y[,k], II_08$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,7), ylim=c(-2.5,2.5))

smooth_data = microRNA(II_05$Y[,k], II_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,7), ylim=c(-2.5,2.5))

smooth_data = microRNA(I_05$Y[,k], I_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,7), ylim=c(-2.5,2.5))

smooth_data = microRNA(I_02$Y[,k], I_02$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,7), ylim=c(-2.5,2.5))

smooth_data = microRNA(Naive$Y[,k], Naive$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-5,7), ylim=c(-2.5,2.5))


######

j = which(mRNA == "hsa.mir.136")
k = which(mRNA == "hsa.mir.376a")
smooth_data = microRNA(II_08$Y[,k], II_08$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-4,4), ylim=c(-2.1,1.8))

smooth_data = microRNA(II_05$Y[,k], II_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-4,4), ylim=c(-2.1,1.8))

smooth_data = microRNA(I_05$Y[,k], I_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-4,4), ylim=c(-2.1,1.8))

smooth_data = microRNA(I_02$Y[,k], I_02$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-4,4), ylim=c(-2.1,1.8))

smooth_data = microRNA(Naive$Y[,k], Naive$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",     
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-4,4), ylim=c(-2.1,1.8))




##################################### Uniqueness by Scenario I

j = which(mRNA == "hsa.mir.148a")
k = which(mRNA == "hsa.mir.210")
smooth_data = microRNA(I_05$Y[,k], I_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2)

######
j = which(mRNA == "hsa.mir.148a")
k = which(mRNA == "hsa.mir.801")

smooth_data = microRNA(I_05$Y[,k], I_05$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,3), ylim=c(-0.5,2.5))

smooth_data = microRNA(I_02$Y[,k], I_02$Y[,j])
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,3), ylim=c(-0.5,2.5))

smooth_data = microRNA(j,k,Sigma_delta = sigma_deltaR1)
plot(smooth_data$x, smooth_data$y,
     type = "l",      
     xlab = paste(mRNA[j]),
     ylab = paste(mRNA[k]), lwd = 2, xlim = c(-3,3), ylim=c(-1.8,1))








ord_1 = order(Z)
smooth_data_1 = lowess(Z[ord_1], pred[ord_1], f = 0.5) 
plot(smooth_data_1$x, smooth_data_1$y,
     type = "l",      
     xlab = "hsa.mir.376a",  #k=14
     ylab = "hsa.mir.377",  #j=15
     lwd = 2#,
#     ylim = c(-1,3.2)
)














