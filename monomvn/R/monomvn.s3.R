#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## summary.monomvn
##
## generic summary method for monomvn class objects,
## records the percentage of marginal (and conditional
## when S=TRUE) columns is printed by analyising S (and Si)

'summary.monomvn' <-
function(object, Si=FALSE, ...)
  {
    rl <- list(obj=object)
    class(rl) <- "summary.monomvn"
    
    ## summary information about the zeros of
    ## the S matrix
    if(!is.null(object$B)) { ## Bayesian version
      
      ## pairwise marginally uncorrelated
      Stri <- 1-object$S.nz[upper.tri(object$S.nz)]
      rl$marg <-  mean(Stri)
      rl$S0 <- nrow(object$S) * apply(1-object$S.nz, 2, mean)
      
      ## conditionally marginally uncorrelated
      if(Si) {
        Sitri <- 1-object$Si.nz[upper.tri(object$Si.nz)]
        rl$cond <- mean(Sitri)
        rl$Si0 <- nrow(object$S) * apply(1-object$Si.nz, 2, mean)
      }
      
    } else { ## MLE/CV version

      ## pairwise marginally uncorrelated
      Stri <- object$S[upper.tri(object$S)]
      rl$marg <-  mean(Stri == 0)
      rl$S0 <- apply(object$S, 2, function(x) { sum(x == 0) })
      
      ## conditionally marginally uncorrelated
      if(Si) {
        Si <- solve(object$S)
        Sitri <- Si[upper.tri(Si)]
        rl$cond <- mean(Sitri == 0)
        rl$Si0 <- apply(Si, 2, function(x) { sum(x == 0) })
      }
    }

    ## print it or return it
    rl
  }


## print.summary.monomvn
##
## print the results of the summary method after first
## calling the print method on the monomvn object.  

'print.summary.monomvn' <-
  function(x, ...)
{
  ## print the monomvn object
  print(x$obj, ...)

  ## count the extra number of things printed
  p <- 0

  ## extra print for the Bayesian case
  if(!is.null(x$obj$B)) cat("On average in the posterior:\n")
  
  ## print the marginally uncorrelated percentage
  if(!is.null(x$marg)) {
    cat(signif(100*x$marg, 3), "% of S is zero", sep="")
    cat(" (pairwise marginally uncorrelated [MUc])\n", sep="")
    if(is.null(x$obj$B)) { ## only in the non-Bayesian case
      m <- sum(x$S0 == length(x$S0) - 1)
      cat("\t", m, " cols (of ",  length(x$S0), " [",
          signif(100*m/length(x$S0), 3),
          "%]) are MI & CI of all others\n", sep="")
    }
    p <- p + 1
  }

  ## print the conditionally uncorrellated percentage
  if(!is.null(x$cond)) {
    cat(signif(100*x$cond, 3), "% of inv(S) is zero", sep="")
    cat(" (pairwise conditionally uncorrelated [CI])\n", sep="")
    p <- p + 1
  }
  
  ## add another newline if we added anything to print.monomvn
  if(p > 0) cat("\n")
}


## plot.summary.monomvn:
##
## make historgrams of the number of zeros in the S matrix
## (and possibly inv(S) matrix) contained in a summary.monomvn
## object 

'plot.summary.monomvn' <-
  function(x, gt0=FALSE, main=NULL, xlab="number of zeros", ...)
{

  ## check if there is anything to plot
  if(all(x$S0 == 0)) {
    cat("S has no zero entries, so there is nothing to plot\n")
    return
  }
  
  ## count the number of things we've plotted
  p <- 0

  ## calculate the dimensions of the plot
  if(!is.null(x$S0) && !is.null(x$Si0))
    par(mfrow=c(1,2))

  ## agument main argument
  smain <- paste(main, "# of zero entries per column", sep="")
  
  ## plot a histogram of the number the zeros for each
  ## asset in S, marking marginal uncorrelation
  if(!is.null(x$S0)) {
    main <- paste(smain, "in S")
    if(gt0) { i <- x$S0 > 0; main <- paste(main, "[>0]") }
    else i <- rep(TRUE, length(x$S0))
    hist(x$S0[i], main=main, xlab=xlab, ...)
    p <- p + 1
  }

  ## plot a histogram of the number of zeros for each
  ## asset in Si, marking conditional uncorrelation
  if(!is.null(x$Si0)) {
    main <- paste(smain, "in inv(S)")
    if(gt0) { i <- x$Si0 > 0; main <- paste(main, "[>0]") }
    else i <- rep(TRUE, length(x$Si0))
    hist(x$Si0[i], main=main, xlab=xlab, ...)
    p <- p + 1
  }
  
  if(p == 0) warning("nothing to plot")
}


## print.monomvn
##
## generic print method for monomvn class objects,
## summarizing the results of a monomvn call

`print.monomvn` <-
function(x, ...)
  {

    ## print information about the call
    cat("\nCall:\n")
    print(x$call)

    ## print information about the methods used
    cat("\nMethods used (p=", x$p, "):\n", sep="")
    um <- sort(unique(x$methods))
    for(u in um) {
      m <- x$methods == u
      cat(sum(m), "\t", u, sep="")
      if(u != "complete" && u != "bcomplete"
         && u != "lsr" && u != "blsr") {
        if(u == "blasso") {
          r <- range(x$lambda2[m])
          ncomp <- "lambda2"
        } else {
          r <- range(x$ncomp[m])
          ncomp <- "ncomp"
        }
        if(u == "ridge") ncomp <- "lambda"
        if(u != "fact") 
          cat(paste(", ", ncomp, " range: [",
                    signif(r[1],5), ",", signif(r[2],5), "]", sep=""))
      }
      cat("\n")
    }
    cat("\n")

    ## in the case of Bayesian regressions
    if(!is.null(x$B)) {
      cat("Bayesian regressions were used with B=", x$B, "\n", sep="")
      cat("burn-in rounds and T=", x$T, " total sampling rounds;\n", sep="")
      cat("see $thin column-wise dynamic thinning level;\n", sep="")
      if(x$rao.s2) cat("Rao-Blackwellized s2 draws were used\n")
      else cat("Standard Park & Casella s2 full-conditional draws were used\n")

      ## say something about DA if any
      if(!is.null(x[['R']])) {
         cat("\nData Augmentation was used for",
             sum(x$R == 2), "entries of y;\n")
         cat("see the $R field for the missingness pattern\n")
      }
      
      ## say something about the reversible jump scheme that was used
      if(x$RJ != "none") {
        cat("\nReversible Jump (RJ) was used to average over\n")
        cat("subsets of columns in the design matrix, for\n")
        if(x$RJ == "p") cat("every parsimonious regression, with a\n")
        else cat("every big-p-small-n regression, with a\n")
        
        ## add in info about the prior
        if(x$mprior[1] == 0)
          cat("uniform prior on m in {1,...,M[i]}\n", sep="")
        else if(length(x$mprior) == 1)
          cat("Bin(m|n=M[i]", ",p=", x$mprior, ") prior\n", sep="")
        else cat("Bin(m|n=M[i],p) prior with p~Beta(",
                 x$mprior[1], ",", x$mprior[2], ")\n", sep="")
      }

      ## say something about the errors used
      if(!is.null(x$theta)) {
        cat("\nStudent-t errors were used with theta=", abs(x$theta), "\n", sep="")
        if(x$theta < 0)
          cat("using a pooled nu for all columns, see $nu\n")
      }
      cat("\n")

      ## check if there are traces
      if(!is.null(x$trace)) {
        cat("Traces are recorded in the $trace field\n\n")
      }

      ## check if there are draws from QP solutions
      if(!is.null(x$QP)) {
        cat("Samples from QP solutions are in the $W field;\n",
            "try summary(obj$W), with monomvn-object \"obj\"\n\n",
            sep="")
      }
    }
  }


## plot.monomvn
##
## generic print method for monomvn class objects,
## summarizing the results of a bmonomvn call
## (only works for bmonomvn at this point)

`plot.monomvn` <-
function(x, which=c("mu", "S", "Snz", "Sinz", "QP"),
         xaxis=c("numna","index"), main=NULL, uselog=FALSE, ...)
{
  ## check that we have bmonomnv (Bayesian) output
  if(is.null(x$B)) stop("plots only supported for bmonomvn output")

  ## check the which argument
  which <- match.arg(which)

  ## check the xaxis argument and construct the x-axis
  xaxis <- match.arg(xaxis)
  if(xaxis == "numna") {
    labs <- x$na[x$o]
    xlab <- "number missing"
  } else {
    labs <- 1:length(x$mu)
    xlab <- "index"
  }

  ## check the uselog argument
  if(length(uselog) > 1 || !is.logical(uselog))
    stop("uselog should be a scalar logical")
  
  ## plot info about mu
  if(which == "mu") {

     ## construct the response
    if(uselog) {
      y <- log(sqrt(x$mu.var[x$o]))
      uselog <- "log"
    } else {
      y <- sqrt(x$mu.var[x$o])
      uselog <- NULL
    }

    ## plot it
    if(is.null(main)) main <- paste(uselog, "sd(mu) by", xlab) 
    plot(labs, y, xlab=xlab, ylab=paste("sd(", which, ")", sep=""))
    title(main)

  } else if(which == "S" || which == "Snz" || which == "Sinz") {
    ## plot info about S

    ## generate unique lables by divding up the unit interval
    ## because image() requires it
    if(xaxis == "numna") {
      d <- labs[duplicated(labs)]
      if(length(d) > 0) {
        for(i in 1:length(d)) {
          rl <- labs == d[i]
          s <- c(0,cumsum(rep(1/(sum(rl)),sum(rl)-1)))
          labs[rl] <- labs[rl] + s
        }
      }
    }

    ## construct the response
    if(which == "S") S <- x$S.var
    else if(which == "Snz") S <- x$S.nz
    else S <- x$Si.nz

    ## use log?
    if(uselog) {
      y <- log(sqrt(S[x$o,x$o]))
      uselog <- "log"
    } else {
      y <- sqrt(S[x$o,x$o])
      uselog <- NULL
    }
    
    ## create main title, image, and plot
    if(is.null(main)) {
      if(which == "S") main <- paste(uselog, "sd(S) by", xlab)
      else if(which == "Snz") main <- paste(uselog, "p(S) != 0 by", xlab)
      else main <- paste(uselog, "p(inv(S)) != 0 by", xlab)
    }
    image(labs, labs, y, xlab=xlab, ylab=xlab)
    title(main)

  } else { ## plot a summary of the weights W

    ## check to make sure QP=TRUE was used to sample from W
    if(is.null(x$W)) stop("must run bmonomvn with QP=TRUE to sample W")

    ## arrange the xaxis differently since the boxplot uses names
    if(xaxis=="numna"){
      labs <- (x$na[-(1:(ncol(x$y)-x$QP$m))])[x$QP$o]
      W <- data.frame(x$W[,x$QP$o])
    } else {
      labs <- ((ncol(x$y)-x$QP$m)+1):ncol(x$y)
      W <- data.frame(x$W)
    }

    ## create a custom x label
    xlab <- paste("ordered by", xlab)
    if(length(labs) > 10) names <- labs
    else names <- names(W)

    ## make the boxplot
    if(is.null(main)) main <- "Boxplots of QP weights w"
    boxplot(W, ylab="w", names=names, xlab=xlab, main=main, ...)

    ## add the MAP point
    points(as.numeric(W[x$which.map,]), col=2, pch=21, cex=1.5)

    ## consider adding a legend, or making logical legend argument
    ## to this function
  }
}
