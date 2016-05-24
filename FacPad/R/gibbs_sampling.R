gibbs_sampling <-
function(matrixY,matrixL,alpha_tau=1,beta_tau=0.01,tau_sig=1,max_iter=10000,thin=10,alpha_sigma=0.7,beta_sigma=0.3,file_name){
  
#library(Rlab)
#library(MASS)


# The input data matrix (G by J, G is the number of genes, J is the number of samples,normalized to mean=0 and sd=1 for each gene)
matrixY<-t(scale(t(matrixY)))
G<-nrow(matrixY)
J<-ncol(matrixY)

# The input binary matrix containing the prior knowledge about the genes associated with each KEGG pathway.
# Each pathway is treated as one latent factor. This matrix has dimension G by K, whereas K is the number of latent factors


K<-ncol(matrixL)



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


matrixW_chain <- matrix(0,nrow=max_iter/thin,ncol=sum(matrixL))
matrixX_chain <- matrix(0,nrow=max_iter/thin,ncol=length(matrixX))
tau_g_chain <- matrix(0,nrow=max_iter/thin,ncol=G)



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
     
      ## Update the values of Wg
         Lg <- matrixL[g,]
         matrixW[g,Lg==0] <- 0
         K1 <- sum(Lg)
         if(K1>1){
            submatrixX <- matrixX[Lg==1,]             
            COV_g <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma * diag(K1)
            }else{
            submatrixX <- t(matrix(matrixX[Lg==1,]))            
            COV_g <-  tau_g[g] * submatrixX %*% t(submatrixX) + tau_sigma
            }
         COV_g_inv <- solve(COV_g)
         mean_g <- tau_g[g]* COV_g_inv %*% submatrixX %*% matrix(matrixY[g,])
         matrixW[g,Lg==1] <- mvrnorm(n = 1, mu=mean_g, Sigma=COV_g_inv)
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
           matrixW_chain[r,] <- matrixW[which(matrixL==1)]
           matrixX_chain[r,] <- as.vector(matrixX)
           tau_g_chain[r,] <- tau_g
           if(tau_sig==0){
              tau_sigma_chain[r] <- tau_sigma
           }
     }
  }
  
  if(tau_sig==0){
    save(matrixW_chain,matrixX_chain,tau_g_chain,tau_sigma_chain,file=file_name)
    }else{
    save(matrixW_chain,matrixX_chain,tau_g_chain,file=file_name)
  }
}

