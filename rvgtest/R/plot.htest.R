##
## Plot results of tests
##
## --------------------------------------------------------------------------

plot.rvgt.htest <- function (x, alpha=0.001, ...)

  ## ------------------------------------------------------------------------
  ## Function for plotting graph of p-values against repetitions
  ## ------------------------------------------------------------------------
  ## x     : object of class "rvgt.htest" containing p-values, or
  ##          list of such objects
  ## alpha : significance level
  ## ------------------------------------------------------------------------
{
  ## check arguments
  if (missing(x) || !is.list(x))
    stop ("Argument 'x' missing or invalid.")

  if (alpha <= 0 || alpha > 0.11)
    stop ("Invalid argument 'alpha'.")

  ## we have two cases:
  ## either result is an object of class "rvgt.htest", or
  ## it is a list of such objects.
  ## The former case is transformed into the latter.
  if (class(x) == "rvgt.htest") {
    res <- list(x)
  }
  else {
    ## here we should check the class of the list members
    if (class(x[[1]]) != "rvgt.htest")
      stop ("Invalid argument 'x'.")
    res <- x
  }
  
  ## -- now creat plot ------------------------------------------------------

  ## parameters for plot
  nplots <- length(res)
  loga <- log10(alpha)

  ## limits for ploting area
  xl <- 0
  yl <- 1.1 * loga
  for (i in 1:nplots) {
    pvals <- res[[i]]$pval
    xl <- max(xl, res[[i]]$r * res[[i]]$n)
    yl <- min(yl, log10(pvals))
  }
  yl <- max(-300,yl)
  
  ## colors
  lc <- rainbow(nplots)
  
  ## create plotting aera with labels
  plot(xl,yl,xlim=c(0,xl),ylim=c(yl,0),type="n",
       xlab="sample size", ylab="log10(pvalue)", ...)

  ## draw lines for each set of p-values
  for (i in 1:nplots) {
    pvals <- res[[i]]$pval
    logp <- pmax(log10(pvals),-300)
    r <- res[[i]]$r
    n <- res[[i]]$n
    tss <- (1:r) * n
    if (length(logp)>1)
      lines(tss,logp,type="l",col=lc[i])
    else
      points(tss[1],logp[1],col=lc[i])
    if (nplots>1)
      text(x=tss[1],y=logp[1],pos=2,col=lc[i],
           label=paste(" [",i,"]",sep=""))
  }
  
  ## line for given significance level
  abline(loga,0,col="red",lwd=2,lty=2)
}

## --------------------------------------------------------------------------

print.rvgt.htest <- function (x, ...) {
  cat("\nrvgtest - test:\n")
  cat("   type           =",x$type,"\n");
  cat("   sample size    =",x$n*x$rep,"\n");
  cat("   # break points =",x$breaks,"\n")
  cat("   p-value        =",x$pval[length(x$pval)],"\n\n")
}

## --------------------------------------------------------------------------
