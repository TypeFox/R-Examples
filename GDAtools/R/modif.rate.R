modif.rate <- function(resmca) {
      Q <- ncol(resmca$call$X)
#      if(attr(resmca,'class')[1]=='multiMCA') Q <- length(resmca$eig[[1]])
      seuil <- 1/Q
      e <- resmca$eig[[1]][resmca$eig[[1]]>=seuil]
      pseudo <- (Q/(Q-1)*(e-seuil))^2
      mrate <- round(pseudo/sum(pseudo)*100,2)
      cum.mrate <- cumsum(mrate)
      return(data.frame(mrate,cum.mrate))
      }
