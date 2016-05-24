gibbs_sampling <-
function(matrixY1,matrixY2,matrixL1,matrixL2,eta0,eta1,alpha_tau=1,beta_tau=0.01,tau_sig=0,max_iter=100000,thin=10,alpha_sigma1=0.7,alpha_sigma2=0.7,
beta_sigma1=0.3,beta_sigma2=0.3,file_name){
  
#library(Rlab)
#library(MASS)


# The input data matrix (J by G, G is the number of genes plus the number of drugs, J is the number of individuals)


matrixY<-list(t(scale(t(matrixY1))),t(scale(t(matrixY2))))
G<-c(nrow(matrixY[[1]]),nrow(matrixY[[2]]))
J<-ncol(matrixY[[1]])

# The input binary matrix containing the prior knowledge about the genes and drugs associated with each KEGG pathway.
# Each pathway is treated as one latent factor. This matrix has dimension G by K, whereas K is the number of latent factors

prior1 <- matrixL1
prior2 <- matrixL2
prior1[matrixL1==0] <- eta0[1]
prior1[matrixL1==1] <- 1-eta1[1]
prior2[matrixL2==0] <- eta0[2]
prior2[matrixL2==1] <- 1-eta1[2]
matrixPi <- list(prior1,prior2)
K<-ncol(matrixPi[[1]])
factor_labels<-c(1:K)

# The G by K loading matrix W and the factor activity matrix X

matrixZ1 <- matrix(0,nrow=G[1], ncol=K)
for(i in 1:G[1]){
   for(j in 1:K){
     matrixZ1[i,j] <- rbern(1,matrixPi[[1]][i,j])
}}

matrixZ2 <- matrix(0,nrow=G[2], ncol=K)
for(i in 1:G[2]){
   for(j in 1:K){
     matrixZ2[i,j] <- rbern(1,matrixPi[[2]][i,j])
}}

matrixZ<-list(matrixZ1,matrixZ2)
rm(matrixZ1)
rm(matrixZ2)

matrixPr1 <- matrix(0,nrow=G[1], ncol=K)
matrixPr2 <- matrix(0,nrow=G[2], ncol=K)
matrixPr<-list(matrixPr1,matrixPr2)
rm(matrixPr1)
rm(matrixPr2)

matrixW1 <- matrix(rnorm(n=G[1]*K, mean = 0, sd = 1), nrow=G[1], ncol=K)
matrixW2 <- matrix(rnorm(n=G[2]*K, mean = 0, sd = 1), nrow=G[2], ncol=K)
matrixW<-list(matrixW1,matrixW2)
rm(matrixW1)
rm(matrixW2)
matrixX <- matrix(rnorm(n=J*K, mean = 0, sd = 1), nrow=K,ncol=J)


tau_g1 <- rep(alpha_tau/beta_tau,G[1])
tau_g2 <- rep(alpha_tau/beta_tau,G[2])
tau_g<-list(tau_g1,tau_g2)
rm(tau_g1)
rm(tau_g2)

if(tau_sig==0){
tau_sigma1 <- alpha_sigma1/beta_sigma1
tau_sigma2 <- alpha_sigma2/beta_sigma2
alpha_sigma<-list(alpha_sigma1,alpha_sigma2)
beta_sigma<-list(beta_sigma1,beta_sigma2)
tau_sigma<-list(tau_sigma1,tau_sigma2)
rm(alpha_sigma1)
rm(alpha_sigma2)
rm(beta_sigma1)
rm(beta_sigma2)
rm(tau_sigma1)
rm(tau_sigma2)
tau_sigma_chain <- list(rep(0,max_iter/thin),rep(0,max_iter/thin))
}else{
tau_sigma<-list(tau_sig,tau_sig)
}



matrixZ_chain <- list(matrix(0,nrow=max_iter/thin,ncol=length(matrixZ[[1]])),matrix(0,nrow=max_iter/thin,ncol=length(matrixZ[[2]])))
matrixW_chain <- list(matrix(0,nrow=max_iter/thin,ncol=length(matrixW[[1]])),matrix(0,nrow=max_iter/thin,ncol=length(matrixW[[2]])))
matrixX_chain <- matrix(0,nrow=max_iter/thin,ncol=length(matrixX))
label_chain <- matrix(0,nrow=max_iter/thin,ncol=K)
tau_g_chain <- list(matrix(0,nrow=max_iter/thin,ncol=G[1]),matrix(0,nrow=max_iter/thin,ncol=G[2]))
matrixPr_chain<-list(matrix(0,nrow=max_iter/thin,ncol=length(matrixZ[[1]])),matrix(0,nrow=max_iter/thin,ncol=length(matrixZ[[2]])))



iter<-0


while(iter < max_iter){


     iter <- iter+1
     print(paste("Working on iteration:",iter))


     if(tau_sig==0){
         ## sample the new precision
         for(r in 1:2){
             n <- sum(matrixZ[[r]])
             new_alpha_sigma <- alpha_sigma[[r]] + n/2
             new_beta_sigma <- beta_sigma[[r]] + sum(matrixW[[r]]^2)/2
             tau_sigma[[r]] <- rgamma(1,shape=new_alpha_sigma, rate=new_beta_sigma)
         }
      }



     for(r in 1:2){

        for(g in 1:G[r]){
     
           for(k in 1:K){  

             #  print(k)       
             Zg <- matrixZ[[r]][g,]
             
             ## Evaluate for Zg,k=1
             Zg[k]<-1
             K1 <- sum(Zg)
             if(K1>1){
                submatrixX <- matrixX[Zg==1,]             
                COV_p1_inv <-  tau_g[[r]][g] * submatrixX %*% t(submatrixX) + tau_sigma[[r]] * diag(K1)
                }else{
                submatrixX <- t(matrix(matrixX[Zg==1,]))            
                COV_p1_inv <-  tau_g[[r]][g] * submatrixX %*% t(submatrixX) + tau_sigma[[r]]
                }
             COV_p1 <- solve(COV_p1_inv)
             mean_p1 <- tau_g[[r]][g]* COV_p1 %*% submatrixX %*% matrix(matrixY[[r]][g,])  
             p_1_exp <- 0.5* t(mean_p1) %*% COV_p1_inv %*% mean_p1

           ## Evaluate for Zg,k=0       
             Zg[k]<-0
             K1 <- sum(Zg)
         #   print(K1)
             if(K1==0){
                ratio_01 <- 0
             }else{
                if(K1>1){
                  submatrixX <- matrixX[Zg==1,]             
                  COV_p0_inv <-  tau_g[[r]][g] * submatrixX %*% t(submatrixX) + tau_sigma[[r]] * diag(K1)
                  }else{               
                  submatrixX <- t(matrix(matrixX[Zg==1,]))            
                  COV_p0_inv <-  tau_g[[r]][g] * submatrixX %*% t(submatrixX) + tau_sigma[[r]]
                  }
                COV_p0 <- solve(COV_p0_inv)
                mean_p0 <- tau_g[[r]][g]* COV_p0 %*% submatrixX %*% matrix(matrixY[[r]][g,])
                p_0_exp <- 0.5* t(mean_p0) %*% COV_p0_inv %*% mean_p0
                ratio_01 <- ((1-matrixPi[[r]][g,k])*sqrt(det(COV_p0)))/(matrixPi[[r]][g,k]*tau_sigma[[r]]*sqrt(det(COV_p1)))*exp(p_0_exp-p_1_exp)
               }
               p_1 <- 1/(1+ratio_01)
               p_0 <- 1-p_1
               matrixPr[[r]][g,k]<-p_1
               matrixZ[[r]][g,k] <- sample(c(0,1),1,prob=c(p_0,p_1))
         }

      ## Update the values of Wg

         Zg <- matrixZ[[r]][g,]
         matrixW[[r]][g,Zg==0] <- 0
         K1 <- sum(Zg)
         if(K1>1){
            submatrixX <- matrixX[Zg==1,]             
            COV_g <-  tau_g[[r]][g] * submatrixX %*% t(submatrixX) + tau_sigma[[r]] * diag(K1)
            }else{
            submatrixX <- t(matrix(matrixX[Zg==1,]))            
            COV_g <-  tau_g[[r]][g] * submatrixX %*% t(submatrixX) + tau_sigma[[r]]
            }
         COV_g_inv <- solve(COV_g)
         mean_g <- tau_g[[r]][g]* COV_g_inv %*% submatrixX %*% matrix(matrixY[[r]][g,])
         matrixW[[r]][g,Zg==1] <- mvrnorm(n = 1, mu=mean_g, Sigma=COV_g_inv)
    #     print(g)

        } 
     }


     COV_x_1_inv <- (t(matrixW[[1]]) %*% diag(tau_g[[1]]) %*% matrixW[[1]] + diag(K))
     COV_x_2_inv <- (t(matrixW[[2]]) %*% diag(tau_g[[2]]) %*% matrixW[[2]] + diag(K))
     COV_x <- solve(COV_x_1_inv + COV_x_2_inv)
     COV_x_1 <- solve(COV_x_1_inv)
     COV_x_2 <- solve(COV_x_2_inv)

     for(j in 1:J){         
         mean_x_1 <- COV_x_1 %*% t(matrixW[[1]]) %*% diag(tau_g[[1]]) %*% matrix(matrixY[[1]][,j])
         mean_x_2 <- COV_x_2 %*% t(matrixW[[2]]) %*% diag(tau_g[[2]]) %*% matrix(matrixY[[2]][,j])
         mean_x <- COV_x %*% (COV_x_1_inv %*% mean_x_1 + COV_x_2_inv %*% mean_x_2)         
         matrixX[,j] <- mvrnorm(n = 1, mu=mean_x, Sigma=COV_x)
     }

     for(r in 1:2){
       SSE<- apply((matrixY[[r]] - matrixW[[r]] %*% matrixX)^2,1,sum)
       tau_g[[r]] <- rgamma(G[r],shape=alpha_tau + J/2, rate=beta_tau + SSE/2)
     }


    ## Permutation step
    for(k in 1:(K-1)){
       for(m in (k+1):K){
          new_old_ratio_1 <- prod(matrixPi[[1]][,m]^(matrixZ[[1]][,k]-matrixZ[[1]][,m]))*prod((1-matrixPi[[1]][,m])^(matrixZ[[1]][,m]-matrixZ[[1]][,k]))*prod(matrixPi[[1]][,k]^(matrixZ[[1]][,m]-matrixZ[[1]][,k]))*prod((1-matrixPi[[1]][,k])^(matrixZ[[1]][,k]-matrixZ[[1]][,m]))
          new_old_ratio_2 <- prod(matrixPi[[2]][,m]^(matrixZ[[2]][,k]-matrixZ[[2]][,m]))*prod((1-matrixPi[[2]][,m])^(matrixZ[[2]][,m]-matrixZ[[2]][,k]))*prod(matrixPi[[2]][,k]^(matrixZ[[2]][,m]-matrixZ[[2]][,k]))*prod((1-matrixPi[[2]][,k])^(matrixZ[[2]][,k]-matrixZ[[2]][,m]))
          new_old_ratio <- new_old_ratio_1*new_old_ratio_2
          prob_old <- 1/(1+new_old_ratio)
          prob_new <- 1-prob_old          
          label <- sample(c(0,1),1,prob=c(prob_new,prob_old))     
          if(label==0){ 
             for(r in 1:2){              
               matrixZ[[r]][,c(k,m)] <- matrixZ[[r]][,c(m,k)]
               matrixW[[r]][,c(k,m)] <- matrixW[[r]][,c(m,k)]
               }
               matrixX[c(k,m),] <- matrixX[c(m,k),]
               factor_labels[c(k,m)]<-factor_labels[c(m,k)]
               print("label switch!")
               }
       }  
    }
 

    r<-as.integer(iter/thin)
    if(r==iter/thin){
        for(t in 1:2){
           matrixZ_chain[[t]][r,] <- as.vector(matrixZ[[t]])
           matrixW_chain[[t]][r,] <- as.vector(matrixW[[t]])
           matrixPr_chain[[t]][r,] <- as.vector(matrixPr[[t]])
           matrixX_chain[r,] <- as.vector(matrixX)
           label_chain[r,]<-factor_labels
           tau_g_chain[[t]][r,] <- tau_g[[t]]
           if(tau_sig==0){
              tau_sigma_chain[[t]][r] <- tau_sigma[[t]]
           }
        }
     }
  }

  
  if(tau_sig==0){
    save(matrixZ_chain,matrixW_chain,matrixX_chain,tau_g_chain,tau_sigma_chain,matrixPr_chain,label_chain,file=file_name)
    }else{
    save(matrixZ_chain,matrixW_chain,matrixX_chain,tau_g_chain,matrixPr_chain,label_chain,file=file_name)
  }
}

