##################################################################
### Examination of Measurement Error Correction for "NP_Graph" ###
##################################################################

library(MASS)

n = 400
p = 50
m = 20
R = 10   # number of repeated measurements
n1 = n
mu = rep(0,p)
Sigma = diag(1,p)

Sigma_eta = diag(0.2,p)
Sigma_delta = diag(0.2,m)

IM1_collect = NULL
Omega1_collect = NULL
IM2_collect = NULL
Omega2_collect = NULL
IM3_collect = NULL
Omega3_collect = NULL
IM4_collect = NULL
Omega4_collect = NULL
IM1 = matrix(0,p,m); IM2 = matrix(0,p,m); IM3 = matrix(0,p,m); IM4 = matrix(0,p,m)
Omega1 = matrix(0,m,m); Omega2 = matrix(0,m,m); Omega3 = matrix(0,m,m); Omega4 = matrix(0,m,m)
N = 20

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



DA1 = NP_Graph(W, Z, sigma_eta = diag(0,p), rho=0.2, sigma_delta = 0,label_name=TRUE)
DA2 = NP_Graph(W, Z, sigma_eta = Sigma_eta, rho=0.2, sigma_delta = 0,label_name=TRUE)
DA3 = NP_Graph(W, Z, sigma_eta = diag(0,p), rho=0.2, sigma_delta = mean(diag(Sigma_delta)),label_name=TRUE)
DA4 = NP_Graph(W, Z, sigma_eta = Sigma_eta, rho=0.2, sigma_delta = mean(diag(Sigma_delta)),label_name=TRUE)


IM1 = IM1 + DA1$importance_matrix
IM2 = IM2 + DA2$importance_matrix
IM3 = IM3 + DA3$importance_matrix
IM4 = IM4 + DA4$importance_matrix

Omega1 = Omega1 + DA1$precision_matrix
Omega2 = Omega2 + DA2$precision_matrix
Omega3 = Omega3 + DA3$precision_matrix
Omega4 = Omega4 + DA4$precision_matrix


IM1_collect[[i]] = DA1$importance_matrix
Omega1_collect[[i]] = DA1$precision_matrix

IM2_collect[[i]] = DA2$importance_matrix
Omega2_collect[[i]] = DA2$precision_matrix

IM3_collect[[i]] = DA3$importance_matrix
Omega3_collect[[i]] = DA3$precision_matrix

IM4_collect[[i]] = DA4$importance_matrix
Omega4_collect[[i]] = DA4$precision_matrix


}

###############################################################
###############################################################

IM1 = IM1/N 
IM1 = matrix(unlist(IM1),p,m)
IM1[which(IM1<2)] = 0
####
Omega1[which(abs(Omega1)<0.2)] = 0

IM2 = IM2/N 
IM2 = matrix(unlist(IM2),p,m)
IM2[which(IM3<2)] = 0
####
Omega2[which(abs(Omega2)<0.2)] = 0

IM3 = IM3/N 
IM3 = matrix(unlist(IM2),p,m)
IM3[which(IM3<2)] = 0
####
Omega3[which(abs(Omega3)<0.2)] = 0


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

#Class(B[1:6,1:6],IM4[1:6,1:6])

######################

#Class(C,Omega4)

######################

var_est = function(IM0, Omega0) {

Measure_B = NULL; Measure_C = NULL 

for(i in 1:N) {
IM = matrix(unlist(IM0[[i]]),p,m)
IM[which(IM<2)] = 0
Omega = matrix(unlist(Omega0[[i]]),m,m)
#Omega[which(Omega<0.1)] = 0
Measure_B = rbind(Measure_B, Class(B,IM)) 
Measure_C = rbind(Measure_C, Class(C,Omega))

}
return(c(apply(Measure_B,2,var)*100, apply(Measure_C,2,var)*100)) 
}

Class(B[1:6,1:6],IM1[1:6,1:6])
Class(C,Omega1)
round(var_est(IM1_collect, Omega1_collect),3)

Class(B[1:6,1:6],IM2[1:6,1:6])
Class(C,Omega2)
round(var_est(IM2_collect, Omega2_collect),3)

Class(B[1:6,1:6],IM3[1:6,1:6])
Class(C,Omega3)
round(var_est(IM3_collect, Omega3_collect),3)

Class(B[1:6,1:6],IM4[1:6,1:6])
Class(C,Omega4)
round(var_est(IM4_collect, Omega4_collect),3)


