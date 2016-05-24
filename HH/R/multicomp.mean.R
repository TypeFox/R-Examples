## calculate multiple comparisons from sufficient statistics
multicomp.mean <- function(group, n, ybar, s, alpha=.05,
                           ylabel="ylabel", focus.name="focus.factor", plot=FALSE,
                           lmat, labels=NULL, ...,
                           df=sum(n) - length(n),
                           sigmahat=(sum((n-1)*s^2) / df)^.5) {

  if.R(r={
    multicomp.default <- NA ## make R-2.6.0dev happy
    stop("multicomp.mean works only in S-Plus.  Use aovSufficient and glht in R.")
  },s={})

  rec.n <- diag(1/n)
  names(ybar) <- group
  if (missing(lmat))
    result <- multicomp.default(ybar,
                                df.residual=df,
                                vmat=rec.n*sigmahat^2,
                                ylabel=ylabel,
                                plot=plot,
                                alpha=alpha,
                                labels=labels,
                                ...)
  else {
    if (missing(labels)) labels <- dimnames(lmat)[[2]]
    result <- multicomp.default(ybar,
                                df.residual=df,
                                vmat=rec.n*sigmahat^2,
                                ylabel=ylabel,
                                plot=plot,
                                alpha=alpha,
                                lmat=lmat,
                                labels=labels,
                                ...)
  }
  oldClass(result) <-  c("multicomp.hh", "multicomp")
  if (!is.null(result$lmat))
    dimnames(result$lmat) <- list(group, dimnames(result$table)[[1]])

  dotdotdot <- list(...)
  if (!is.null(dotdotdot$crit.point) & !is.null(dotdotdot$method))
    result$method <- dotdotdot$method
  
  result$focus.name <- focus.name
  result
}
## trace(multicomp.mean, exit=browser)



multicomp.mmc.mean <- function(group, n, ybar, s, ylabel, focus.name,
                               lmat,
                               ...,
                               comparisons="mca",
                               lmat.rows=seq(length=length(ybar)),
                               ry,
                               plot=TRUE,
                               crit.point,
                               iso.name=TRUE,
                               estimate.sign=1,
                               x.offset=0,
                               order.contrasts=TRUE,
                               method="tukey",
                               df=sum(n)-length(n),
                               sigmahat=(sum((n-1)*s^2)/df)^.5) {
  
  if.R(r=stop("multicomp.mmc.mean works only in S-Plus.  Use aovSufficient and mmc in R."),
       s={})

  ## pairwise differences
  if (missing(crit.point)) {
    mca <- multicomp.mean(group, n, ybar, s, ylabel=ylabel, focus.name=focus.name,
                          comparisons=comparisons, plot=FALSE, method=method,
                          df=df, sigmahat=sigmahat, ...)
    crit.point <- mca$crit.point
  }
  else
    mca <- multicomp.mean(group, n, ybar, s, ylabel=ylabel, focus.name=focus.name,
                          comparisons=comparisons,
                          crit.point=crit.point, plot=FALSE, method=method,
                          df=df, sigmahat=sigmahat, ...)
  if (estimate.sign != 0) mca <- multicomp.reverse(mca, estimate.sign)
  
  ## group means
  none <- multicomp.mean(group, n, ybar, s, ylabel=ylabel, focus.name=focus.name,
                         comparisons="none",
                         crit.point=mca$crit.point, plot=FALSE, method=method,
                         df=df, sigmahat=sigmahat, ...)
  nte <- none$table[,"estimate"]
  none$height <- nte * 2
  mca$height <- (nte %*% abs(mca$lmat[lmat.rows,]))[1,]
  
  ## user-specified lmat
  if (!missing(lmat)) {
    lmat <- sweep(lmat, 2,
                     apply(abs(lmat[lmat.rows, , drop=FALSE]), 2, sum)/2, "/")
    lmat.multicomp <- multicomp.mean(group, n, ybar, s, ylabel=ylabel,
                           comparisons="none",
                           crit.point=mca$crit.point,
                           lmat=lmat, plot=FALSE, method=method,
                           df=df, sigmahat=sigmahat, ...)
    if (!is.null(lmat.multicomp$message)) stop(lmat.multicomp$message)
    if (estimate.sign != 0) lmat.multicomp <- multicomp.reverse(lmat.multicomp, estimate.sign)
    lmat.multicomp$height <- (nte %*% abs(lmat.multicomp$lmat))[1,]
  }
  
  if (order.contrasts) {
    mca <- multicomp.order(mca)
    none <- multicomp.order(none)
    if (!missing(lmat)) lmat.multicomp <- multicomp.order(lmat.multicomp)
  }
  
  ## result
  result <- list(mca=mca, none=none)
  if (!missing(lmat)) result$lmat.multicomp <- lmat.multicomp
  if (!missing(ry)) result$ry <- ry
  oldClass(result) <- c("mmc.multicomp", "list")
  
  ## plot
  if (plot) plot.mmc.multicomp(result, iso.name=iso.name, x.offset=x.offset)
  
  return(result)
}

## trace(multicomp.mmc.mean, exit=browser)
## source(hh("splus.library/multicomp.mmc.mean.s"))
