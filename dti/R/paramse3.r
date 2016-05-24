gethseqfullse3 <- function (kstar, gradstats, kappa=NULL, vext = c(1, 1)) 
{
  ngrad <- dim(gradstats$bghat)[2]
  h  <- vr  <- matrix(0,ngrad,kstar)
  if(is.null(kappa)) kappa <- pi/ngrad
  #    if(length(kappa)<kstar) kappa <- rep(kappa[1],kstar)
  prt0 <- Sys.time()
  cat("get sequences of bw, kappa up to kstar=", kstar, " ")
  n <- 0
  for(i in 1:ngrad){
    z <- .Fortran("ghfse3i",
                  as.integer(i),#i4
                  as.integer(kstar),#kstar
                  as.double(gradstats$k456),
                  as.integer(ngrad),
                  as.double(kappa),#kappa
                  as.double(vext),#vext
                  h=double(kstar),
                  vr=double(kstar),#
                  n=integer(1),#
                  as.integer(gradstats$dist),
                  PACKAGE="dti")[c("h","vr","n")]
    h[i,] <- z$h
    vr[i,] <- z$vr 
    n <- n+z$n
    cat(".")
  }
  cat("\n number of positive weights:",n,"mean maximal bandwidth",signif(mean(h[,kstar]),3),"time elapsed:", 
      format(difftime(Sys.time(), prt0), digits = 3), 
      "\n")
  list(h=h,kappa=kappa,vred=vr,n=n)
}
reduceparam <- function(param){
  ind <- param$ind[4,]==param$ind[5,]
  param$ind <- param$ind[,ind]
  param$w <- param$w[ind]
  param$n <- sum(ind)
  h <- max(param$h)
  rind <- param$ind[1,]*h*2*h*2+param$ind[2,]*h*2+param$ind[3,]
  oind <- order(rind)
  param$ind <- param$ind[,oind]
  param$w <- param$w[oind]
  starts <- cumsum(rle(rind[oind])$lengths)
  param$nstarts <- length(starts)
  param$starts <- c(0,starts)
  param
}


