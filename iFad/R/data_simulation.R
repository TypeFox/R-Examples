data_simulation <-
function(K,G1,G2,J,eta0,eta1,density,alpha_tau=1,beta_tau=0.01,SNR=0,file_name){

#  library(Rlab)
#  library(MASS)

  matrixL1<-matrix(0,G1,K)
  matrixL2<-matrix(0,G2,K)

  n1<-round(length(matrixL1)*density[1])
  n2<-round(length(matrixL2)*density[2])

  matrixL1[sample(length(matrixL1),n1)]<-1
  matrixL2[sample(length(matrixL2),n2)]<-1

  # The G by K Bernoulli prior matrix Pie1 and Pie2
  matrixPi1 <- matrix(0,nrow=G1,ncol=K)
  matrixPi2 <- matrix(0,nrow=G2,ncol=K)
  matrixPi1[matrixL1==0] <- eta0[1]
  matrixPi1[matrixL1==1] <- 1-eta1[1]
  matrixPi2[matrixL2==0] <- eta0[2]
  matrixPi2[matrixL2==1] <- 1-eta1[2]

  matrixZ1 <- matrix(0,nrow=G1, ncol=K)
  for(i in 1:G1){
    for(j in 1:K){
     matrixZ1[i,j] <- rbern(1,matrixPi1[i,j])
  }}

  matrixZ2 <- matrix(0,nrow=G2, ncol=K)
  for(i in 1:G2){
    for(j in 1:K){
      matrixZ2[i,j] <- rbern(1,matrixPi2[i,j])
  }}


  matrixW1 <- matrix(rnorm(n=G1*K, mean = 0, sd = 1), nrow=G1, ncol=K)
  matrixW1[matrixZ1==0]<-0
  matrixW2 <- matrix(rnorm(n=G2*K, mean = 0, sd = 1), nrow=G2, ncol=K)
  matrixW2[matrixZ2==0]<-0
  matrixX<-matrix(rnorm(n=J*K, mean = 0, sd = 1), nrow=K,ncol=J)
  Y1_mean <- matrixW1 %*% matrixX
  Y2_mean <- matrixW2 %*% matrixX

  if(SNR==0){
      tau_g1 <- rgamma(G1,shape=alpha_tau, rate=beta_tau)
      tau_g2 <- rgamma(G2,shape=alpha_tau, rate=beta_tau)
      sigma1 <- diag(1/tau_g1)
      sigma2 <- diag(1/tau_g2)
    }else{
      sigma1 <- diag(rep(K/SNR,G1))
      sigma2 <- diag(rep(K/SNR,G2))
    }  

  matrixY1 <- matrix(0,nrow=G1,ncol=J)  
  matrixY2 <- matrix(0,nrow=G2,ncol=J)

  for(i in 1:J){
    matrixY1[,i] <- mvrnorm(n = 1, mu = Y1_mean[,i], Sigma = sigma1, tol = 1e-6, empirical = FALSE)
    matrixY2[,i] <- mvrnorm(n = 1, mu = Y2_mean[,i], Sigma = sigma2, tol = 1e-6, empirical = FALSE)
  }

  save(matrixL1,matrixL2,matrixPi1,matrixPi2,matrixW1,matrixW2,matrixX,matrixY1,matrixY2,matrixZ1,matrixZ2,sigma1,sigma2,Y1_mean,Y2_mean,file=file_name) 

}

