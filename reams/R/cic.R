cic <-
function(y, X, nperms=499, covests=NULL, nullcic=NULL) {
    require(leaps)
    n = length(y)
    sigma2.hat = crossprod(lsfit(X, y)$resid) / -diff(dim(X))
    sigma2.y = var(y)
    if (is.null(nullcic)) nullcic=as.numeric((n-1)*sigma2.y+2*sigma2.hat)/n
    realleaps = leaps(X, y, names=names(X), method='r2')
    if (is.null(covests))  {
        permleaps <- vector("list", nperms)
        sumcovs = matrix(NA, nperms, ncol(X))
        for (pp in 1:nperms) {
  	        if (pp/20 == floor(pp/20)) cat('Permutation', pp, '\n')
  	        perm <- sample(n)
  	        permleaps[[pp]] = leaps(X, y[perm], names=names(X), method='r2', nbest=1)
  	        for (k in 1:ncol(X)) {
  	  	        X.k = X[ , permleaps[[pp]]$which[k, ]]
  	  	        eta.k = y[perm] - lsfit(X.k, y[perm])$resid
  	  	        sumcovs[pp, k] = crossprod(scale(y[perm], scale=FALSE), eta.k)
  	        }
        }
        covests = colMeans(sumcovs)
    }
    else permleaps = NULL
    enp = covests/sigma2.y+1
    covests.rep = rep(covests, table(realleaps$size))
    cic1 = (n-1) * sigma2.y * (1-realleaps$r2) / n
    cic2 = 2 * sigma2.hat * covests.rep / (n * sigma2.y)
    cic3 = 2 * sigma2.hat / n
    cic = cic1 + cic2 + cic3
    list(realleaps=realleaps, permleaps=permleaps, covests=covests, enp=enp, cic=cic, nullcic=nullcic, best=realleaps$which[which.min(cic), ])
}

