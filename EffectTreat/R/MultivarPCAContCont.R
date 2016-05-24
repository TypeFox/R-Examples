Multivar.PCA.ContCont <- function(Sigma_TT, Sigma_TS, Sigma_SS, T0T1=seq(-1, 1, by=.01)) { 
  
  Results <- T0T1_results <- NULL
  
  for (i in 1: length(T0T1)) {   

    Sigma_TT_hier <- Sigma_TT
    Sigma_TT_hier[1,2] <- Sigma_TT_hier[2, 1] <- (sqrt(Sigma_TT[1,1])*sqrt(Sigma_TT[2,2])) * T0T1[i]
    
    Sigma_temp_a <- 
      rbind(cbind(Sigma_TT_hier, Sigma_TS), 
      cbind(t(Sigma_TS), Sigma_SS))
    Cor_c <- cov2cor(Sigma_temp_a)
    a_1 <- matrix(c(-1, 1), nrow=1)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)  
    
    if (Min.Eigen.Cor > 0) {
      
      PCA <- (a_1 %*% Sigma_TS %*% solve(Sigma_SS) %*% t(Sigma_TS) %*% t(a_1))/
                (a_1 %*% Sigma_TT_hier %*% t(a_1))
                 
    Results <- rbind(Results, PCA)
    T0T1_results <- rbind(T0T1_results, T0T1[i])
    }
  }
  
  if (is.null(Results)==FALSE){  
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  colnames(Results) <- "PCA"
  Pos.Def <- nrow(Results)
  }
  
  if (is.null(Results)==TRUE){  
    cat("No solutions found, try specifying another range of rho_T0T1 values. \n")
    rownames(Results) <- NULL
    Pos.Def <- 0
    T0T1_results <- NA
  }
  
  fit <- 
    list(Total.Num.Matrices=length(T0T1), Pos.Def=Pos.Def, T0T1=as.numeric(t(T0T1_results)), 
         PCA=Results$PCA, Call=match.call())
  
  class(fit) <- "Multivar.PCA.ContCont"
  fit
}