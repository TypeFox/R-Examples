Pos.Def.Matrices <- function(T0T1=seq(0, 1, by=.2), T0S0=seq(0, 1, by=.2), T0S1=seq(0, 1, by=.2), 
                        T1S0=seq(0, 1, by=.2), T1S1=seq(0, 1, by=.2), S0S1=seq(0, 1, by=.2)){

combins <- expand.grid(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1)

Generated.Matrices <- matrix(NA, 1, 8)
colnames(Generated.Matrices) <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "Min.Eigen.Sigma", "Pos.Def.Status")
Generated.Matrices <- na.exclude(Generated.Matrices)

for (i in 1: nrow(combins)) {
 T0T1 <- combins[i, 1]  
 T0S0 <- combins[i, 2]  
 T0S1 <- combins[i, 3]
 T1S0 <- combins[i, 4]
 T1S1 <- combins[i, 5]  
 S0S1 <- combins[i, 6]  
    
 Sigma_c <- diag(4)   
 Sigma_c[lower.tri(Sigma_c, diag=FALSE)] <- matrix(c(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1))   
 Sigma_c[upper.tri(Sigma_c, diag=FALSE)] <- matrix(c(T0T1, T0S0, T1S0, T0S1, T1S1, S0S1))   
 Min.Eigen.Sigma <- try(min(eigen(Sigma_c)$values), TRUE)   

 if (Min.Eigen.Sigma > 0) {
  PD_status <- c(1) 
  results <- cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1, Min.Eigen.Sigma, PD_status)
  Generated.Matrices <- rbind(Generated.Matrices, results)        
    }
 if (Min.Eigen.Sigma <= 0) {
   PD_status <- c(0) 
   results <- cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1, Min.Eigen.Sigma, PD_status)
   Generated.Matrices <- rbind(Generated.Matrices, results)        
 }
}

Generated.Matrices <<- data.frame(Generated.Matrices)

fit <- list(Generated.Matrices=Generated.Matrices)

}