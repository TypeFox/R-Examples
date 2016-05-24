sp.assoc <-
function(matr)
{  
   if(any(is.na(matr))){
       matr <- na.omit(matr)
       cat("The data matrix contains NA, and have been removed.\n")
   }
   matr[matr>1] <- 1
   frq <- freq.calc(matr)
   N <- nrow(matr) ##Number of plots
   S <- ncol(matr) ##Number of species
   n <- apply(matr, 2, sum) ## Number of plots occupied by certain species
   Tj <- apply(matr, 1, sum) ## Total number of species for each plot
   t <- mean(Tj) ## Mean species number for all the plots
   STsq <- (1/N)*sum((Tj-t)^2) ## Variance of species number
   sigmaTsq <- sum((1-frq)*frq) ## Variance of species relative frequency
   VR <- STsq/sigmaTsq ##Variance ratio: VR>1 Positively associated, VR<1 Negative associated
   W <- VR * N  # W statistic value, used in comparison with Chi(n).
   result <- list(pi = frq, N = N, S = S, Tj = Tj, Numspmean = t,
   sigmaTsq = sigmaTsq, STsq = STsq, var.ratio = VR, W = W)
   return(result)
}

