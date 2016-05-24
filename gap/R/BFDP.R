BFDP <- function(a,b,pi1,W,logscale=FALSE)
{
  if (logscale) {
     thetahat <- a
     V <- b
     PO <- (1-pi1)/pi1
  } else {
     thetahat <- log(a)
     V <- ((log(b)-log(a))/1.96)^2
     pi0 <- 1-pi1
     PO <- pi0/(1-pi0)
  }
  postvar <- V + W
  r <- W / postvar
  z <- thetahat / sqrt(V)
  ABF <- exp(-r*z^2/2) / sqrt(1-r)
  ABFDP <- ABF*PO/(ABF*PO+1)
  pH0 <- dnorm(thetahat,mean=0,sd=sqrt(V))
  pH1 <- dnorm(thetahat,mean=0,sd=sqrt(postvar))
  BF <- pH0/pH1
  BFDP <- BF*PO/(BF*PO+1)
  invisible(list(pH0=pH0,pH1=pH1,BF=BF,BFDP=BFDP,ABF=ABF,ABFDP=ABFDP))
}

# 31-7-2007
