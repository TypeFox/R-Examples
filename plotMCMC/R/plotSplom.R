plotSplom <- function(mcmc, axes=FALSE, between=0, div=1, log=FALSE, base=10, ...)
{
  ## 1  Parse args
  if(is.mcmc.list(mcmc))
    mcmc <- as.mcmc(mcmc)
  if(is.mcmc(mcmc))
    mcmc <- as.data.frame(mcmc)
  ellipsis <- as.list(substitute(list(...)))[-1]
  if(is.null(dim(mcmc)))
    stop("Argument 'mcmc' must contain more than one chain, arranged in columns.")

  ## 2  Transform
  mcmc <- if(log) log(mcmc/div,base=base) else mcmc/div

  ## 3  Draw plot
  if(!axes && is.null(ellipsis$oma))
    suppressWarnings(pairs(mcmc, gap=between, oma=c(0,0,0,0), xaxt="n", yaxt="n", ...))
  else
    suppressWarnings(pairs(mcmc, gap=between, ...))

  invisible(NULL)
}
