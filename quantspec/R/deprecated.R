################################################################################
# This code is from version 0.2; don't use these functions any more!
# Version 1.0-0 was completely redesigned and offers the same functionality.
#
# The old codes is soley kept for compatability reasons and to warn users about
# its depreciation. It will be removed in the near future.
################################################################################

################################################################################
# This code used to be in ct.R
################################################################################

#' Defunct functions in package quantspec
#'
#' These functions have been declared defunct since Version 1.0-1.
#'
#' @name quantspec-defunct
#' @aliases ct
#' @keywords Defunct
#'
#' @param i1 Parameter of DEFUNCT function.
#' @param i2 Parameter of DEFUNCT function.
#' @param n   Parameter of DEFUNCT function.
#' @param X Parameter of DEFUNCT function.
#' @param taus Parameter of DEFUNCT function.
#' @param omegas Parameter of DEFUNCT function.
#' @param fromRanks Parameter of DEFUNCT function.
#' @param showProgressBar Parameter of DEFUNCT function.
#' @param LPG Parameter of DEFUNCT function.
#' @param F Parameter of DEFUNCT function.
#' @param CL Parameter of DEFUNCT function.
#' @param hRange Parameter of DEFUNCT function.
#' @param hOffset Parameter of DEFUNCT function.
#' @param ylabel Parameter of DEFUNCT function.
#' @param oma Parameter of DEFUNCT function.
#' @param mar Parameter of DEFUNCT function.
#' @param cex.lab Parameter of DEFUNCT function.
#' @param W Parameter of DEFUNCT function.

ct <- function (i1, i2, n) {
  .Defunct("getValues")

#  i1*(n+0.5)-n+i2-0.5*i1^2
}

################################################################################
# This code used to be in LaplacePeriodogram.R
################################################################################

#' @name quantspec-defunct
#' @aliases LaplacePeriodogram
LaplacePeriodogram <-
    function (X, taus, omegas=1:(ceiling(length(X)/2)-1), fromRanks=TRUE, showProgressBar=FALSE) {

  .Defunct("quantilePG")

#  # Verify if all parameters are valid
#  if (!is.vector(X) & !is.ts(X) & !is.zoo(X)) {
#    stop("'X' needs to be specified as a vector, a ts or a zoo object")
#  }
#
#  if (is.ts(X)) {
#    X <- X[1:(length(X))]
#  }
#
#  if (is.zoo(X)) {
#    X <- coredata(X)
#  }
#
#  if (!(is.vector(taus) && is.numeric(taus) && prod(taus>0) && prod(taus<1))) {
#    stop("'taus' needs to be specified as a vector of quantile orders")
#  }
#
#  is.wholenumber <-
#      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
#
#  if (!(is.vector(omegas) && prod(is.wholenumber(omegas)))) {
#    stop("'omegas' needs to be specified as a vector of integers")
#  } else {
#    if (!(prod(omegas>=0) * prod(omegas<=floor(length(X)/2)))) {
#      stop(paste("all values in omegas need to be integeres between 1 and",floor(length(X)/2)))
#    }
#  }
#
#  # Convert Data to standardized ranks
#  if (fromRanks) {
#    X <- rank(X)/length(X)
#  }
#
#  # Initialize Progress Bar
#  if (showProgressBar) {
#    maxPb <- length(omegas)+1
#    #pb <- tkProgressBar(title = "Calculate Laplace Periodogram", min=0, max=maxPb, width = 400)
#  }
#
#  # Define and initiate auxiliary variables.
#  n <- length(X)
#  ltau <- length(taus)
#  res <- rep(floor(n/2)*ltau*(ltau+1)/2,0)
#  ctau <- 1
#  ires <- 1
#
#  # Then calculate the Laplace Periodogram at all frequencies in omegas.
#  for (i in omegas) {
#    if (showProgressBar) {
#      #iPb <- getTkProgressBar(pb)+1
#      #setTkProgressBar(pb, iPb, label=paste(round(100*iPb/maxPb,0),"% done"))
#    }
#    # Define the harmonic regressors.
#    omega <- (2 * pi * i) / n
#    D <- cos(omega*1:n)
#    S <- sin(omega*1:n)
#
#    # Then performe the quantile regression.
#    qregSol <- coef(rq(X ~ 1 + D + S, taus))
#    if (length(taus)>1) {
#      qregSol <- qregSol[2:3,]
#    } else {
#      qregSol <- matrix(qregSol[2:3],ncol=1)
#    }
#
#    # Then for each combination of taus ...
#    for (i1 in 1:ltau) {
#      for (i2 in i1:ltau) {
#        B1 <- qregSol[,i1]
#        B2 <- qregSol[,i2]
#        res[ires] <- n/4 * complex(real = t(B1) %*% B2, imaginary = B1[2]*B2[1]-B1[1]*B2[2])
#        ires <- ires+1
#      }
#    }
#  }
#
#  if (showProgressBar) {
#    #close(pb)
#  }
#
#  # At last prepare the row and column labels
#  cnames <- rep(ltau*(ltau+1)/2,0)
#  for (i1 in 1:ltau) {
#    for (i2 in i1:ltau) {
#      cnames[ctau] <- paste(taus[i1],taus[i2],sep="-")
#      ctau <- ctau+1
#    }
#  }
#  rnames=2*pi/n*omegas
#
#  # Then return the result as a matrix.
#  matrix(res,ncol=(ltau*(ltau+1)/2),nrow=length(omegas),byrow=T,dimnames=list(rnames,cnames))
}

################################################################################
# This code used to be in plotLaplacePeriodogram.R
################################################################################

#' @name quantspec-defunct
#' @aliases plotLaplacePeriodogram
# @importFrom rje is.subset
plotLaplacePeriodogram <-
    function(LPG, taus, F=1:length(LPG[,1]), CL=1:length(taus),
        hRange=FALSE, hOffset=FALSE,
        ylabel=expression({{hat(f)}[n]^{list(tau[1],tau[2])}}(omega)),
        oma=c(2.5,2.5,2.5,2.5),
        mar=c(4.5,4.5,1,0)+0.1,
        cex.lab=1.5) {

  .Defunct("plot")
  
#  warning("This function has been depreciated since Version 1.0\n",
#      "It will be removed in the near future!")
#
#  # Verify if all parameters are valid
#
#  if (!(is.vector(taus) && is.numeric(taus) && prod(taus>0) && prod(taus<1))) {
#    stop("'taus' needs to be specified as a vector of quantile orders")
#  }
#
#  if (!(dim(LPG)[2] == (length(taus)+1)*length(taus)/2)) {
#    stop("length of vector 'taus' does not correspond to number of columns in LPG")
#  }
#
#  is.wholenumber <-
#      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
#
#  if (!(is.vector(CL) && prod(is.wholenumber(CL)) && prod(CL > 0) && prod(CL <= length(taus)))) {
#    stop(paste("'CL' needs to be specified as a vector of indices between 1 and",length(taus)))
#  }
#
#  if (!(is.vector(F) && is.integer(F))) {
#    stop("'F' needs to be specified as a vector of integers")
#  } else {
#    omegas <- round(as.numeric(rownames(LPG))*dim(LPG)[1]/pi)
#    if (!is.subset(F,omegas)) {
#      stop("'F' needs to be a subset of 'omegas' in 'LPG'")
#    }
#  }
#
#  # Define and initiate auxiliary variable.
#  ltau <- length(taus)
#
#  # Set graphics parameters.
#  par(mfcol=c(length(CL),length(CL)),oma=oma,mar=mar,cex.lab=cex.lab)
#
#  # Create Matrix of values to plot
#  # Used to determine xlim and ylim
#  M <- matrix(c(Re(LPG[F,]),Im(LPG[F,])),nrow=length(F))
#
#  # For each combination i1, i2 create a plot
#  for (i1 in CL) {
#    for (i2 in CL) {
#
#      # Labels only for the plots at left and bottom
#      if (i1 == CL[1]) {
#        yl <- substitute(paste(tau[1],"=",k),list(k=taus[i2]))
#      } else {yl <- ""}
#      if (i2 == CL[length(CL)]) {
#        xl <- substitute(paste(tau[2],"=",k),list(k=taus[i1]))
#      } else {xl <- ""}
#
#      # Arrange points to be ploted as X,Y matrix
#      # Frequencies to be plotted
#      if (i2 >= i1) {
#        ctau <- ct(i1,i2,ltau)
#        # Real parts of the estimate
#        Est <- matrix(c(as.numeric(rownames(LPG[F,]))/(2*pi),Re(LPG[F,ctau])),ncol=2)
#      } else {
#        ctau <- ct(i2,i1,ltau)
#        # Imaginary parts of the estimate
#        Est <- matrix(c(as.numeric(rownames(LPG[F,]))/(2*pi),Im(LPG[F,ctau])),ncol=2)
#        ctau <- ctau+length(LPG[1,])
#      }
#
#      # Determine ylimits
#      if (hOffset) {
#        ylimits <- c(min(M),max(M))
#      } else {
#        if (hRange) {
#          range <- max(apply(M,2,max)-apply(M,2,min))-max(M[,ctau])+min(M[,ctau])
#          ylimits <- c(min(M[,ctau])-range/2,max(M[,ctau])+range/2)
#        } else {
#          ylimits <- c(min(M[,ctau]),max(M[,ctau]))
#        }
#      }
#
#      xmin <- as.numeric(rownames(LPG)[min(F)])/(2*pi)
#      xmax <- as.numeric(rownames(LPG)[max(F)])/(2*pi)
#      plot(x=c(xmin,xmax),y=ylimits,type="n", xlab=xl, ylab=yl)
#      lines(Est)
#    }
#  }
#  mtext(expression(omega),outer=TRUE,side=1,line=0)
#  mtext(ylabel,outer=TRUE,side=2,line=0)
#
#  par(mfcol=c(1,1),oma=oma,mar=c(5, 4, 4, 2) + 0.1,cex.lab=1)
}

################################################################################
# This code used to be in smoothedLaplacePeriodogram.R
################################################################################

#' @name quantspec-defunct
#' @aliases smoothedLaplacePeriodogram
smoothedLaplacePeriodogram <-
    function(LPG, taus, W) {

  .Defunct("smoothedPG")

#  # Verify if all parameters are valid
#
#  if (!is.vector(W) && (length(W) %% 2) == 0) {
#    if (is.tskernel(W)) {
#      W <- W$coef
#      W <- c(W[length(W):2],W[1:length(W)])
#    } else {
#      stop("'W' needs to be specified as a vector with an odd number of entries or as a tskernel object")
#    }
#  }
#
#  if (!(is.vector(taus) && is.numeric(taus) && prod(taus>0) && prod(taus<1))) {
#    stop("'taus' needs to be specified as a vector of quantile orders")
#  }
#
#  if (!(dim(LPG)[2] == (length(taus)+1)*length(taus)/2)) {
#    stop("length of vector 'taus' does not correspond to number of columns in LPG")
#  }
#
#  # Define and initiate auxiliary variables.
#  lomega <- length(LPG[,1])
#  ltau <- length(taus)
#  lW <- length(W)
#  res <- matrix(rep(0,(lomega*ltau*(ltau+1)/2)),ncol=(ltau*(ltau+1)/2))
#  # Normalize vector of weights
#  W <- W / sum(W)
#
#  # Extend LPG in the beginning and at the end using the symmetry.
#  m <- (lW-1)/2
#  extendedLPG <- matrix(c(t(Conj(LPG[(m-1):1,])),Conj(LPG[1,]),t(LPG[1:lomega,]),
#          t(Conj(LPG[(lomega):(lomega-m),]))),ncol=(ltau*(ltau+1)/2),byrow=T)
#
#  # Calculate the weighted averages for each frequency.
#  for (i in 1:lomega) {
#    res[i,] <- W %*% extendedLPG[i:(i+2*m),]
#  }
#
#  # Add column and row names and return the result.
#  colnames(res) <- colnames(LPG)
#  rownames(res) <- rownames(LPG)
#  res
}
