Sim.Data.Counterfactuals <- function(N.Total=2000, mu_c=c(0, 0, 0, 0), T0S0=0, T1S1=0, T0T1=0, T0S1=0, T1S0=0, S0S1=0, Seed=sample(1:1000, size=1)) {  

if ((length(T0T1)>1 | length(T0S1)>1 |  length(T1S0)>1 | length(S0S1)>1)==TRUE) stop("Please make sure that T0T1, T0S1, T1S0, and S0S1 are scalar.")
cors_nok <- (T1S0 < -1 | T1S0 > 1 | T0S1 < -1 | T0S1 > 1 | T0T1 < -1 | T0T1 > 1 | S0S1 < -1 | S0S1 > 1) 
if (cors_nok==TRUE) stop("Please make sure that T0T1, T0S1, T1S0, and S0S1 are within the [-1; 1] range.")

Sigma_c <- diag(4)        
Sigma_c[lower.tri(Sigma_c, diag=FALSE)] <- matrix(c(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1))  
Sigma_c[upper.tri(Sigma_c, diag=FALSE)] <- matrix(c(T0T1, T0S0, T1S0, T0S1, T1S1, S0S1))  
Min.Eigen.Sigma <- try(min(eigen(Sigma_c)$values), TRUE)    

if (Min.Eigen.Sigma <= 0) stop("The minimum eigenvalue of matrix Sigma_c is not positive definite (i.e., Sigma_c is not a valid correlation matrix). Please specify other values for T0S1, T1S0, T1S1, and/or S0S1.")

try(U <- chol(Sigma_c), TRUE)   
set.seed(Seed)
random <- replicate(4, runif(N.Total))    
Ideal <- random%*%U    
Z <- rbinom(N.Total,1,0.5)  
Z[Z==0] <- c(-1)
Ideal <- cbind(Ideal, Z)
colnames(Ideal) <- c("T0", "T1", "S0", "S1", "Z")
Ideal <- data.frame(Ideal)
mu_c <- mu_c-.5
Ideal$S0 <- Ideal$S0 + mu_c[1]
Ideal$S1 <- Ideal$S1 + mu_c[2]
Ideal$T0 <- Ideal$T0 + mu_c[3] 
Ideal$T1 <- Ideal$T1 + mu_c[4]
Data.Counterfactuals <- NULL
Data.Counterfactuals <<- Ideal 

fit <- list(Data.Counterfactuals=Data.Counterfactuals)

}
