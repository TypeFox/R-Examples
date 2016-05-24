### tuning parameter selection for function zipath
### output: return of the best model 
tuning.zipath <- function(formula, data, weights, subset, na.action, offset, standardize=TRUE,
                   family = c("poisson", "negbin", "geometric"),penalty = c("enet", "mnet", "snet"), 
lambdaCountRatio = .0001, lambdaZeroRatio = c(.1, .01, .001), maxit.theta=1, gamma.count=3, gamma.zero=3, ...){
    fam <- match.arg(family)
    pen <- match.arg(penalty)
    ratio <- lambdaZeroRatio
        fm.ratio <- vector("list", length=length(ratio))
        reskk <- rep(NA, length=length(ratio))
        for(kk in 1:length(ratio)){
            fm.ratio[[kk]] <- zipath(formula, data = data, family = fam, gamma.count=gamma.count, gamma.zero=gamma.zero, nlambda=100, lambda.count.min.ratio=lambdaCountRatio, lambda.count=NULL, lambda.zero= NULL, lambda.zero.min.ratio=ratio[kk], maxit.theta=maxit.theta, theta.fixed=FALSE, trace=FALSE, penalty=pen, rescale=FALSE)
            reskk[kk] <- min(fm.ratio[[kk]]$bic)
        }
      fm <- fm.ratio[[which.min(reskk)]]
    return(fm)
}
