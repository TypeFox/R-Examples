PCA.ContCont <- function(T0S, T1S, T0T0=1, T1T1=1, SS=1, T0T1=seq(-1, 1, by=.01)) { 
  
  Results <- na.exclude(matrix(NA, 1, 6))  
  colnames(Results) <- c("T0T1", "T0S", "T1S", "PCA", "Sigma.Delta.T", "delta") 
  combins <- expand.grid(T0T1, T0S, T1S)
  
  for (i in 1: nrow(combins)) {   
    T0T1 <- combins[i, 1]   #corr
    T0S <- combins[i, 2]   #corr
    T1S <- combins[i, 3]   #corr

    Sigma_c <- diag(3)         
    Sigma_c[2,1] <- Sigma_c[1,2] <- T0T1 * (sqrt(T0T0)*sqrt(T1T1))
    Sigma_c[3,1] <- Sigma_c[1,3] <- T0S * (sqrt(T0T0)*sqrt(SS))
    Sigma_c[3,2] <- Sigma_c[2,3] <- T1S * (sqrt(T1T1)*sqrt(SS))
    Sigma_c[1,1] <- T0T0
    Sigma_c[2,2] <- T1T1
    Sigma_c[3,3] <- SS
    Cor_c <- cov2cor(Sigma_c)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)  
    
    if (Min.Eigen.Cor > 0) {
      PCA <- 
        ((sqrt(T1T1)*T1S) - (sqrt(T0T0)*T0S))/
        sqrt(T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * T0T1))
           
      if ((is.finite(PCA))==TRUE){
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * T0T1)
        delta <- (1-(PCA**2)) * sigma.delta.T
        results.part <- as.vector(cbind(T0T1, T0S, T1S, PCA, sigma.delta.T, delta))
        Results <- rbind(Results, results.part)
        rownames(Results) <- NULL}
    }
  }
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  Total.Num.Matrices <- nrow(combins)
  
  fit <- 
    list(Total.Num.Matrices=Total.Num.Matrices, Pos.Def=Results[,1:3], PCA=Results$PCA, GoodSurr=Results[,4:6], Call=match.call())
  
  class(fit) <- "PCA.ContCont"
  fit
  
}