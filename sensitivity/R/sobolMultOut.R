# Aggregated Sobol' indices for multidimensional outputs
# Author: Bertrand Iooss, 2015

# Refs:
# M. Lamboni,H. Monod and D. Makowski, Multivariate sensitivity analysis to measure global contribution of input factors in dynamic models,
#         Reliability Engineering and System Safety, 96:450-459, 2011.
# F. Gamboa, A. Janon, T. Klein and A. Lagnoux, Sensitivity indices for multivariate outputs,
#         Electronic Journal of Statistics, 8:575-603, 2014

sobolMultOut <- function(model = NULL, q = 1, X1, X2, MCmethod = "sobol", 
                         plotFct=FALSE, ...) 
{
# q = dimension of the output vector
# Mcmethod: limited to the following Sobol' estimators: 
#   sobol, sobol2002, sobol2007, soboljansen, sobolmartinez, sobolGP, sobolmara
# PlotFct: plot of the Sobol' indices for each output (for 1D functional output)

  # WARNINGS : no bootstrap ; no tell function ; do not work with sobolEff (pb with the output format $V)
  
#  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) 
#        stop("The samples X1 and X2 must have the same dimensions")
  p <- ncol(X1)

  Tot=FALSE
  if(MCmethod=="sobol2002"||MCmethod=="sobol2007"||MCmethod=="soboljansen"
     ||MCmethod=="sobolmartinez"||MCmethod=="sobolGP") Tot=TRUE

  if(MCmethod == "sobol") x0 <- sobol(model=NULL, X1, X2, ...)
  if(MCmethod == "sobol2002") x0 <- sobol2002(model=NULL, X1, X2, ...)
  if(MCmethod == "sobol2007") x0 <- sobol2007(model=NULL, X1, X2, ...)
  if(MCmethod == "soboljansen") x0 <- soboljansen(model=NULL, X1, X2, ...)
  if(MCmethod == "sobolmartinez") x0 <- sobolmartinez(model=NULL, X1, X2, ...)
#  if(MCmethod == "sobolEff") x0 <- sobolEff(model=NULL, X1, X2, ...)
  if(MCmethod == "sobolGP") x0 <- sobolGP(model=NULL, X1, X2, ...)  
  if(MCmethod == "sobolmara") x0 <- sobolmara(model=NULL, X1, ...)  
  
  x <- x0
  xx <- as.list(c()) # liste concatenant tous les objets Sobol pour toutes les sorties
  
  if (!is.null(model)){
    y <- model(x0$X, ...)

    V <- 0
    for (i in 1:q){
      tell(x0,y[,i])
      V <- V + x0$V[1]
      if (plotFct==TRUE) xx[[i]] <- x0
    }

    x <- x0
    x$V[1] <- V
    for (i in 1:p){
      x$S[i,1] <- V[1+i,1] / V[1,1]
      if (Tot == T) x$T[i,1] <- V[p+1+i,1] / V[1,1]
    }

    if (plotFct == TRUE){
      S <- matrix(0,nrow=q,ncol=p)
      if (Tot == T) ST <- matrix(0,nrow=q,ncol=p)
      for (i in 1:q){
        S[i,] <- xx[[i]]$S[1:p,]
        if (Tot == T) ST[i,] <- xx[[i]]$T[1:p,]
      }
      
      if (Tot == T) par(mfrow=c(2,1))
      plot(0,ylim=c(0,1),xlim=c(1,q),main="First order Sobol indices",ylab="")
      for (i in 1:p) lines(S[,i],col=i)
      legend(x = "topright", legend = dimnames(X1)[[2]], lty=1, col=1:p, cex=0.6)
      if (Tot == T){
        plot(0,ylim=c(0,1),xlim=c(1,q),main="Total Sobol indices",ylab="")
        for (i in 1:p) lines(ST[,i],col=i)
        legend(x = "topright", legend = dimnames(X1)[[2]], lty=1, col=1:p, cex=0.6)
      }
    }
  }
  
  return(x)

}

# TODO : tell function


print.sobolMultOut <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nFirst order indices:\n")
    print(x$S)
    cat("\nTotal indices:\n")
    print(x$T)
  }
}


plot.sobolMultOut <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:p)+.3, add = TRUE)
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}
