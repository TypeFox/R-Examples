plotCumu <- function(mcmc, probs=c(0.025,0.975), div=1, log=FALSE, base=10, main=NULL, xlab="Iterations", ylab="Value",
                     lty.median=1, lwd.median=2, col.median="black", lty.outer=2, lwd.outer=1, col.outer="black", ...)
{
  ## 1  Parse args
  ellipsis <- as.list(substitute(list(...)))[-1]
  probs <- c(min(probs), 0.5, max(probs))
  lty <- c(lty.outer, lty.median)
  lwd <- c(lwd.outer, lwd.median)
  col <- c(col.outer, col.median)

  ## 2  Transform
  mcmc <- if(log) log(mcmc/div,base=base) else mcmc/div

  ## 3  Draw plot
  if(is.null(dim(mcmc)) && !is.list(mcmc))  # vector/mcmc
  {
    if(is.null(ellipsis$ann) && is.null(ellipsis$axes))    # -,-
      cumuplot(mcmc(mcmc), probs=probs, ask=FALSE, lty=lty, lwd=lwd, col=col, ann=FALSE, axes=FALSE, ...)
    if(is.null(ellipsis$ann) && !is.null(ellipsis$axes))   # -,axes
      cumuplot(mcmc(mcmc), probs=probs, ask=FALSE, lty=lty, lwd=lwd, col=col, ann=FALSE, ...)
    if(!is.null(ellipsis$ann) && is.null(ellipsis$axes))  # ann,-
      cumuplot(mcmc(mcmc), probs=probs, ask=FALSE, lty=lty, lwd=lwd, col=col, axes=FALSE, ...)
    if(!is.null(ellipsis$ann) && !is.null(ellipsis$axes))  # ann,axes
      cumuplot(mcmc(mcmc), probs=probs, ask=FALSE, lty=lty, lwd=lwd, col=col, ...)
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
    cumuplot(mcmc(mcmc), probs=probs, ask=FALSE, xlab=xlab, ylab=ylab, ...)
  }

  invisible(NULL)
}
