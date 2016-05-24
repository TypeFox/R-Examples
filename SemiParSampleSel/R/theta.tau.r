theta.tau <- function( BivD, theta.star ) {

  if(BivD %in% c("N", "AMH") )    { theta <- tanh(theta.star); 	   theta1 <- abs(mean(theta)); KendTau <- tau(normalCopula(theta1))  }
  if(BivD %in% c("C0", "C180") )  { theta <- exp(theta.star);	   theta1 <- abs(mean(theta)); KendTau <- tau(claytonCopula(theta1)) }
  if(BivD %in% c("C90","C270") )  { theta <- -exp(theta.star);	   theta1 <- abs(mean(theta)); KendTau <- -tau(claytonCopula(theta1))}
  if(BivD %in% c("J0", "J180"))   { theta <- 1+exp(theta.star);	   theta1 <- abs(mean(theta)); KendTau <- tau(joeCopula(theta1))     }
  if(BivD %in% c("J90", "J270"))  { theta <- -(1+exp(theta.star)); theta1 <- abs(mean(theta)); KendTau <- -tau(joeCopula(theta1))    }
  if(BivD=="FGM")                 { theta <- tanh(theta.star);	   theta1 <- abs(mean(theta)); KendTau <- tau(fgmCopula(theta1))     }
  if(BivD=="F")                   { theta <- theta.star;           theta1 <- abs(mean(theta)); KendTau <- tau(frankCopula(theta1))   }
  if(BivD=="AMH")                 { theta <- tanh(theta.star);	   theta1 <- abs(mean(theta)); KendTau <- tau(amhCopula(theta1))     }
  if(BivD %in% c("G0", "G180"))   { theta <- 1+exp(theta.star);	   theta1 <- abs(mean(theta)); KendTau <- tau(gumbelCopula(theta1))  }
  if(BivD %in% c("G90", "G270"))  { theta <- -(1+exp(theta.star)); theta1 <- abs(mean(theta)); KendTau <- -tau(gumbelCopula(theta1)) }
  
  #KendTau <- NULL
  
  cbind( theta , KendTau)

} 