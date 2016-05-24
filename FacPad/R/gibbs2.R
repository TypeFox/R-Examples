gibbs2 <-
function(matrixY,matrixL,eta0,eta1,alpha_tau=1,beta_tau=0.01,tau_sig=1,max_iter=10000,thin=10,alpha_sigma=0.7,beta_sigma=0.3,file_name){
  
#library(Rlab)
#library(MASS)


# The input data matrix (G by J, G is the number of genes, J is the number of samples,normalized to mean=0 and sd=1 for each gene)
matrixY<-t(scale(t(matrixY)))
G<-nrow(matrixY)
J<-ncol(matrixY)

# The input binary matrix containing the prior knowledge about the genes associated with each KEGG pathway.
# Each pathway is treated as one latent factor. This matrix has dimension G by K, whereas K is the number of latent factors


K<-ncol(matrixL)
matrixPi <- matrixL
matrixPi[matrixL==0] <- eta0
matrixPi[matrixL==1] <- 1-eta1

matrixZ <- matrix(0,nrow=G, ncol=K)
for(i in 1:G){
   for(j in 1:K){
     matrixZ[i,j] <- rbern(1,matrixPi[i,j])
}}

matrixPr <- matrix(0,nrow=G, ncol=K)

# The G by K loading matrix W and the factor activity matrix X

matrixW <- matrix(rnorm(n=G*K, mean = 0, sd = 1), nrow=G, ncol=K)

matrixX <- matrix(rnorm(n=J*K, mean = 0, sd = 1), nrow=K,ncol=J)


tau_g <- rep(alpha_tau/beta_tau,G)


if(tau_sig==0){
tau_sigma <- alpha_sigma/beta_sigma
tau_sigma_chain <- rep(0,max_iter/thin)
}else{
tau_sigma<-tau_sig
}


matrixW_chain <- matrix(0,nrow=max_iter/thin,ncol=length(matrixL))
matrixZ_chain <- matrix(0,nrow=max_iter/thin,ncol=length(matrixL))
matrixX_chain <- matrix(0,nrow=max_iter/thin,ncol=length(matrixX))
tau_g_chain <- matrix(0,nrow=max_iter/thin,ncol=G)
matrixPr_chain<-matrix(0,nrow=max_iter/thin,ncol=length(matrixL))



iter<-0


while(iter < max_iter){


     iter <- iter+1
     print(paste("Working on iteration:",iter))


     if(tau_sig==0){
         ## sample the new precision

             n <- sum(matrixL)
             new_alpha_sigma <- alpha_sigma + n/2
             new_beta_sigma <- beta_sigma + sum(matrixW^2)/2
             tau_sigma <- rgamma(1,shape=new_alpha_sigma, rate=new_beta_sigma)
         
      }
         


     for(g in 1:G){

         for(k in 1:K){ 
           if( matrixPi[g,k]>0 & matrixPi[g,k]<1 ){     
             Zg <- matrixZ[g,] 
             K1 <- sum(Zg) 
             if( K1<= 1 ){
                ratio_01 <- 0
             }else{
                ## Evaluate for Zg,k=0       
                Zg[k]<-0
                K1 <- sum(Zg) 
                if(K1>1){
                  submatrixX <- matrixX[Zg==1,]             
                  COV_p0_inv <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma * diag(K1)
                  }else{               
                  submatrixX <- t(matrix(matrixX[Zg==1,]))            
                  COV_p0_inv <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma
                  }
                COV_p0 <- solve(COV_p0_inv)
                mean_p0 <- tau_g[g]* COV_p0 %*% submatrixX %*% matrix(matrixY[g,])
                p_0_exp <- 0.5* t(mean_p0) %*% COV_p0_inv %*% mean_p0

                ## Evaluate for Zg,k=1
                Zg[k]<-1
                K1 <- sum(Zg)
                if(K1>1){
                  submatrixX <- matrixX[Zg==1,]             
                  COV_p1_inv <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma * diag(K1)
                  }else{
                  submatrixX <- t(matrix(matrixX[Zg==1,]))            
                  COV_p1_inv <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma
                  }
                COV_p1 <- solve(COV_p1_inv)
                mean_p1 <- tau_g[g]* COV_p1 %*% submatrixX %*% matrix(matrixY[g,])  
                p_1_exp <- 0.5* t(mean_p1) %*% COV_p1_inv %*% mean_p1

                ratio_01 <- ((1-matrixPi[g,k])*sqrt(det(COV_p0)))/(matrixPi[g,k]*tau_sigma*sqrt(det(COV_p1)))*exp(p_0_exp-p_1_exp)
               }

          p_1 <- 1/(1+ratio_01)
          p_0 <- 1-p_1
          matrixPr[g,k]<-p_1
          matrixZ[g,k] <- sample(c(0,1),1,prob=c(p_0,p_1))

          }else{     
             matrixPr[g,k] <- matrixPi[g,k]
             matrixZ[g,k] <- matrixPi[g,k]
             }
         }

         ## Update the values of Wg
         Zg <- matrixZ[g,]
         matrixW[g,Zg==0] <- 0
         K1 <- sum(Zg)
         if(K1>1){
            submatrixX <- matrixX[Zg==1,]             
            COV_g <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma * diag(K1)
            }else{
            submatrixX <- t(matrix(matrixX[Zg==1,]))            
            COV_g <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma
            }
         COV_g_inv <- solve(COV_g)
         mean_g <- tau_g[g]* COV_g_inv %*% submatrixX %*% matrix(matrixY[g,])
         matrixW[g,Zg==1] <- mvrnorm(n = 1, mu=mean_g, Sigma=COV_g_inv)
#        print(paste("Working on gene:",g))

     } 

     COV_x_inv <- (t(matrixW) %*% diag(tau_g) %*% matrixW + diag(K))
     COV_x <- solve(COV_x_inv)

     mean_x <- COV_x %*% t(matrixW) %*% diag(tau_g) %*% matrixY
     mvgen<-function(vector){
        t<-mvrnorm(n = 1, mu=vector, Sigma=COV_x)
        return(t)
     }
     matrixX<-apply(mean_x,2,mvgen)

     SSE<- apply((matrixY - matrixW %*% matrixX)^2,1,sum)
     tau_g <- rgamma(G,shape=alpha_tau + J/2, rate=beta_tau + SSE/2)
     
 

     r<-as.integer(iter/thin)
     if(r==iter/thin){
           matrixW_chain[r,] <- as.vector(matrixW)
           matrixZ_chain[r,] <- as.vector(matrixZ)
           matrixPr_chain[r,] <- as.vector(matrixPr)
           matrixX_chain[r,] <- as.vector(matrixX)
           tau_g_chain[r,] <- tau_g
           if(tau_sig==0){
              tau_sigma_chain[r] <- tau_sigma
           }
     }
  }
  
  if(tau_sig==0){
    save(matrixW_chain,matrixZ_chain,matrixPr_chain,matrixX_chain,tau_g_chain,tau_sigma_chain,file=file_name)
    }else{
    save(matrixW_chain,matrixZ_chain,matrixPr_chain,matrixX_chain,tau_g_chain,file=file_name)
  }
}




















