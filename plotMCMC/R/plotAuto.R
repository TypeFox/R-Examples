plotAuto <- function(mcmc, thin=1, log=FALSE, base=10, main=NULL, xlab="Lag", ylab="Autocorrelation", lty=1, lwd=1,
                     col="black", ...)
{
  ## 1  Parse args
  ellipsis <- as.list(substitute(list(...)))[-1]

  ## 2  Transform
  mcmc <- if(log) log(mcmc,base=base) else mcmc

  ## 3  Draw plot
  if(is.null(dim(mcmc)) && !is.list(mcmc))  # vector/mcmc
  {
    if(is.null(ellipsis$ann) && is.null(ellipsis$axes))    # -,-
      autocorr.plot(window(mcmc(mcmc),thin=thin), lty=lty, lwd=lwd, col=col, ask=FALSE, ann=FALSE, axes=FALSE, ...)
    if(is.null(ellipsis$ann) && !is.null(ellipsis$axes))   # -,axes
      autocorr.plot(window(mcmc(mcmc),thin=thin), lty=lty, lwd=lwd, col=col, ask=FALSE, ann=FALSE, ...)
    if(!is.null(ellipsis$ann) && is.null(ellipsis$axes))   # ann,-
      autocorr.plot(window(mcmc(mcmc),thin=thin), lty=lty, lwd=lwd, col=col, ask=FALSE, axes=FALSE, ...)
    if(!is.null(ellipsis$ann) && !is.null(ellipsis$axes))  # ann,axes
      autocorr.plot(window(mcmc(mcmc),thin=thin), lty=lty, lwd=lwd, col=col, ask=FALSE, ...)
    if(is.null(ellipsis$ann) || ellipsis$ann)
    {
      title(main=main, ...)
      title(xlab=xlab, ...)
      title(ylab=ylab, ...)
    }
    if(is.null(ellipsis$axes) || ellipsis$axes)
    {
      axis(1, ...)
      axis(2, ...)
      box()
    }
  }
  else  # matrix/data frame/list/mcmc.list
  {
    mcmc <- if(is.mcmc.list(mcmc)) as.data.frame(sapply(mcmc,unclass)) else as.data.frame(mcmc)
    autocorr.plot(window(mcmc(mcmc),thin=thin), lty=lty, lwd=lwd, col=col, ask=FALSE, ...)
  }

  invisible(NULL)
}
