plot.BEST <-
function(x, which=c("mean", "sd", "effect", "nu"), credMass=0.95,
                    ROPE=NULL, compVal=0, showCurve=FALSE, ...) {
  # This function plots the posterior distribution for one selected item. 
  # Description of arguments:
  # x is mcmc.list object of the type returned by function BESTmcmc.
  # which indicates which item should be displayed; possible values are "mean", "sd",
  #  "effect" or "nu".
  # ROPE is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE.
  # compVal is a scalar specifying the value for comparison.
  # showCurve if TRUE the posterior should be displayed as a fitted density curve
  #   instead of a histogram (default).

  # TODO additional sanity checks.
  # Sanity checks:
  if(!inherits(x, "data.frame"))
    stop("x is not a valid BEST object")
  if(ncol(x) == 3 && all(colnames(x) == c("mu","nu","sigma"))) {
    oneGrp <- TRUE
  } else if (ncol(x) == 5 && all(colnames(x) == c("mu1", "mu2","nu","sigma1","sigma2"))) {
    oneGrp <- FALSE
  } else {
    stop("x is not a valid BEST object")
  }

  # Deal with ... argument
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]

  whichID <- match.arg(which)

  toPlot <- switch(whichID,
    mean = if(oneGrp) x$mu else x$mu1 - x$mu2,
    sd = if(oneGrp) x$sigma else x$sigma1 - x$sigma2,
    effect = if(oneGrp) (x$mu - compVal) / x$sigma else
       (x$mu1 - x$mu2) /
          sqrt( ( x$sigma1^2 + x$sigma2^2 ) / 2 ),
    nu = log10(x$nu) )

  if(is.null(dots$main))
    dots$main <- switch(whichID,
      mean = if(oneGrp) "Mean" else "Difference of Means",
      sd = if(oneGrp) "Standard Deviation" else "Difference of Std. Dev.s",
      effect = "Effect Size",
      nu = "Normality")

  if(is.null(dots$xlab))
    dots$xlab <- switch(whichID,
      mean = if(oneGrp) bquote(mu) else bquote(mu[1] - mu[2]),
      sd = if(oneGrp) bquote(sigma) else bquote(sigma[1] - sigma[2]),
      effect = if(oneGrp) bquote( (mu-.(compVal)) / sigma ) else
        bquote( (mu[1]-mu[2]) / sqrt((sigma[1]^2 +sigma[2]^2 )/2 ) ),
      nu = bquote("log10("*nu*")"))

  if(whichID=="nu" && !is.null(compVal) && compVal == 0)
    compVal <- NULL
  if(whichID=="sd" && oneGrp && !is.null(compVal) && compVal == 0)
    compVal <- NULL
  # Plot posterior distribution of selected item:
  histinfo <- plotPost(toPlot, credMass=credMass, ROPE=ROPE, showCurve=showCurve,
                  showMode=whichID != "mean",
                  compVal=compVal, graphPars=dots) 

  return(invisible(histinfo))
}
