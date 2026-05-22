#################################################
### Constants off-diag enties in ME variances ###
#################################################
library(MASS)
n = 400
p = 50
m = 20
R = 10   # number of repeated measurements
n1 = n
mu = rep(0,p)
Sigma = diag(1,p)

Sigma_eta = matrix(0.2*0.8,p,p) + diag(0.2*0.2,p)
Sigma_delta = matrix(0.2*0.8,m,m) + diag(0.2*0.2,m)

IM4_collect = NULL
Omega4_collect = NULL
IM4 = matrix(0,p,m)
Omega4 = matrix(0,m,m)

N = 100

for(i in 1:N) {


X = mvrnorm(n,mu,Sigma)

Y1 = sin(X[,1]) + rnorm(n,0,1)
Y2 = sin(X[,1]) + cos(X[,2]) + sin(X[,3]) + rnorm(n,0,1)
Y3 = exp(Y1) + X[,1] + X[,2]^2 + X[,4]^2 + rnorm(n,0,1)
Y4 = sin(X[,1]) + cos(X[,4]) + rnorm(n,0,1)
Y5 =  rnorm(n,0,1)
Y6 = cos(Y4) +Y5 + exp(X[,3]) + X[,4]^2 + rnorm(n,0,1)
Y = mvrnorm(n,rep(0,(m-6)),diag(1,(m-6)))
Y = cbind(Y1,Y2,Y3,Y4,Y5,Y6,Y)

################
#B = matrix(0,p,m)
#B[1,1]=1; B[c(1:3),2]=1; B[c(1,2,4),3]=1; B[c(1,4),4]=1; B[c(3,4),6]=1

#C = matrix(0,m,m)
#C[1,3]=1; C[3,1]=1; C[c(4,5),6]=1; C[6,c(4,5)]=1
#C = C + diag(2,m)
#Y = X%*%B + mvrnorm(n,rep(0,m),solve(C))
#################

W = Y + mvrnorm(n,rep(0,m),Sigma_delta)
Z = X + mvrnorm(n,rep(0,p),Sigma_eta)


#### repeated measurements ####
Xval = mvrnorm(n1,mu,Sigma)

Y1 = sin(Xval[,1]) + rnorm(n1,0,1)
Y2 = sin(Xval[,1]) + cos(Xval[,2]) + sin(Xval[,3]) + rnorm(n1,0,1)
Y3 = exp(Y1) + Xval[,1] + Xval[,2]^2 + Xval[,4]^2 + rnorm(n1,0,1)
Y4 = sin(Xval[,1]) + cos(Xval[,4]) + rnorm(n1,0,1)
Y5 =  rnorm(n1,0,1)
Y6 = cos(Y4) +Y5 + exp(Xval[,3]) + Xval[,4]^2 + rnorm(n1,0,1)
Yval = mvrnorm(n1,rep(0,(m-6)),diag(1,(m-6)))
Yval = cbind(Y1,Y2,Y3,Y4,Y5,Y6,Yval)
#########################
#Yval = Xval%*%B + mvrnorm(n1,rep(0,m),solve(C))
#########################
Wsum = matrix(0,n1,m); Zsum =  matrix(0,n1,p); Wval = NULL; Zval = NULL
for(r in 1:R) {

Wval[[r]] = Yval + mvrnorm(n1,rep(0,m),Sigma_delta)
Zval[[r]] = Xval + mvrnorm(n1,rep(0,p),Sigma_eta)
Wsum = Wsum + Wval[[r]]
Zsum = Zsum + Zval[[r]]

}
Sigma_delta_hat = matrix(0,m,m); Sigma_eta_hat = matrix(0,p,p)
for(r in 1:R) {
Sigma_delta_hat = Sigma_delta_hat + t(Wval[[r]] - Wsum/R) %*% ((Wval[[r]] - Wsum/R))
Sigma_eta_hat = Sigma_eta_hat + t(Zval[[r]] - Zsum/R) %*% ((Zval[[r]] - Zsum/R))
}


Sigma_delta_hat = Sigma_delta_hat/(R*(n1-1))
Sigma_eta_hat = Sigma_eta_hat/(R*(n1-1))

DA4 = NP_Graph(W, Z, sigma_eta = Sigma_eta_hat, rho=0.2, sigma_delta = mean(diag(Sigma_delta_hat)),
label_name = TRUE)



IM4 = IM4 + DA4$importance_matrix
Omega4 = Omega4 + DA4$precision_matrix
IM4_collect[[i]] = DA4$importance_matrix
Omega4_collect[[i]] = DA4$precision_matrix


}

###############################################################
###############################################################

IM4 = IM4/N 
IM4 = matrix(unlist(IM4),p,m)
IM4[which(IM4<2)] = 0
####
Omega4[which(abs(Omega4)<0.2)] = 0


###########################################

B = matrix(0,p,m)
B[1,1]=1; B[c(1:3),2]=1; B[c(1,2,4),3]=1; B[c(1,4),4]=1; B[c(3,4),6]=1

C = matrix(0,m,m)
C[1,3]=1; C[3,1]=1; C[c(4,5),6]=1; C[6,c(4,5)]=1
C = C + diag(2,m)


Class = function(Act,Pred) {
TP = 0; FP = 0; TN = 0; FN = 0
p = dim(Act)[1]; m = dim(Act)[2]
for(i in 1:p) {
for(j in 1:m) {
TP = TP + (Act[i,j]!=0 & Pred[i,j] !=0) *1
FP = FP + (Act[i,j]==0 & Pred[i,j] !=0) *1
TN = TN + (Act[i,j]==0 & Pred[i,j] ==0) *1
FN = FN + (Act[i,j]!=0 & Pred[i,j] ==0) *1
}
}
Spe = TN / (TN+FP)
Sen = TP / (TP+FN)
Mcc = (TN*TP - FN*FP) / sqrt((TP+FN) * (TN+FP) * (TP+FP) * (TN+FN) )

return(round(c(Spe,Sen,Mcc),3))
}

Class(B[1:6,1:6],IM4[1:6,1:6])

######################

Class(C,Omega4)

######################

Measure_B = NULL; Measure_C = NULL 

for(i in 1:N) {
IM = matrix(unlist(IM4_collect[[i]]),p,m)
IM[which(IM<2)] = 0
Omega = matrix(unlist(Omega4_collect[[i]]),m,m)
#Omega[which(Omega<0.1)] = 0
Measure_B = rbind(Measure_B, Class(B,IM)) 
Measure_C = rbind(Measure_C, Class(C,Omega))

}


Class(B[1:6,1:6],IM4[1:6,1:6])
Class(C,Omega4)
apply(Measure_B,2,var)*100
apply(Measure_C,2,var)*100

