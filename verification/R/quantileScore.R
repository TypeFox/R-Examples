`quantileScore` <-
function (obs, pred, p, breaks, ...) 
{
    id <- is.finite(obs) & is.finite(pred)
    obs <- obs[id]
    pred <- pred[id]
    #pred <- round(pred, 8)
    #breaks <- round(breaks, 8 )

     # baseline
    obar <- quantile(obs,p,type=8)

     # overall quantile score original forecasts
    qs.orig <- mean( check.func(obs-pred,p) )

     # discretize forecast values
    XX <- quantile2disc(pred, bins = breaks)
    pred <- XX$new
    y.i <- XX$mids

     # number of forecast-observation pairs
    N <- length(obs)
    K <- length(y.i)

     # number of forecasts within each bin
    N.pred <- aggregate(pred, by = list(pred), length)$x

     # conditional observed quantile
    obar.i <- aggregate(obs, by = list(pred), quantile, p)$x

     # overall quantile score
    qs <- mean( check.func(obs-pred,p) )
    qs.baseline <- mean( check.func(obs-obar,p) )
    ss <- 1 - qs/qs.baseline

     # decomposition
    d.CQ <- rep(NA,K)
    d.PQ <- rep(NA,K)
    for(k in 1:K){
      ind.k <- which(pred==as.factor(y.i[k]))
      d.CQ[k] <- sum(check.func(obs[ind.k]-obar,p) - check.func(obs[ind.k]-obar.i[k],p))
      d.PQ[k] <- sum(check.func(obs[ind.k]-y.i[k],p) - check.func(obs[ind.k]-obar.i[k],p))
    }
  
    qs.rel <- sum( d.PQ )/N
    qs.res <- sum( d.CQ )/N
    qs.uncert <- qs.baseline

    check <- qs-(qs.rel - qs.res + qs.uncert)
    prob.y <- N.pred/N

    out <- list(qs.orig=qs.orig, qs = qs, qs.baseline = qs.baseline,
        ss = ss, qs.reliability = qs.rel, qs.resol = qs.res, 
        qs.uncert = qs.uncert, y.i = y.i, obar.i = obar.i, prob.y = prob.y, 
        obar = obar, breaks = breaks, check = check)

    class(out) <- "quantile"

    return(out)
}

quantile2disc <- function(x, bins) {
  ## converts continuous (quantile) forecasts into a range of discrete
  ## (quantile) forecasts assigned to the mean value within each bin
 
  if(max(x) >  max(bins) | min(x) < min(bins) ){stop("
  Bins must span the interval of predictions.") }

  xx  <- cut(x, breaks = bins, include.lowest = TRUE)
  ind <- seq(1,length(bins)-1)[xx]
  
  mids <- aggregate(x, by = list(ind), mean)
  
  #falls ein bin nicht besetzt ist
  new.mids <- rep(NA,length(bins)-1)
  new.mids[mids$Group.1] <- mids$x

  return(list(new = new.mids[xx], mids = mids$x))
}

check.func <- function(u,p) {
  ## calculates check function for values u and given quantile level p
  ## Yu et al. (2001)

  rho <- (abs(u) + (2*p - 1)*u)*0.5

}
