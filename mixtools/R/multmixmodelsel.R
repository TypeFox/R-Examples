multmixmodel.sel <- function (y, comps = NULL, ...) 
{
    if (class(y)=="list" && !is.null(y$y)) {
      y <- y$y
    }
    n = dim(y)[1]
    p = dim(y)[2]
    m = min(apply(y, 1, sum))
#    m = unique(apply(y, 1, sum))
#    if (length(m) > 1) {
#        stop("Each row of y must have same total number of observations")
#    }
    max.allowed.comp = floor((m + 1)/2)
    if (is.null(comps)) 
        comps = 1:max.allowed.comp
    if (max(comps) > max.allowed.comp) {
        stop(paste("No more than", max.allowed.comp, "components allowed",
                   "with", m, "multinomial trials"))
    }
    aic = NULL
    bic = NULL
    caic = NULL
    icl = NULL
    ll = NULL
    theta = matrix(0, 0, p)
    lambda = NULL
    for (k in sort(comps)) {
#        cat("Testing", k, "components:  ")
#        newrows = k - nrow(theta)
        tmp <- multmix.init(y, k = k)
        theta <- tmp$theta
        lambda <- tmp$lambda
      if (k!=1){
      em = multmixEM(y, lambda = lambda, theta = theta, k = k, ...)
        loglik = em$loglik
        lambda = em$lambda
        theta = em$theta
#        cat(em$iter, "iterations.\n")
      } else loglik = sum(log(exp(apply(y,1,ldmult,theta=theta))))
        aic = c(aic, loglik - (p * k - 1))
        bic = c(bic, loglik - log(n) * (p * k - 1)/2)
    caic = c(caic, loglik - (log(n) + 1) * (p * k - 1)/2)
    if (k==1) {
    icl = c(icl, loglik - log(n) * (p * k - 1)/2)
    } else icl = c(icl, loglik - log(n) * (p * k - 1)/2 - sum(lambda * log(lambda)))
        ll = c(ll, loglik)
    }
    out = rbind(aic, bic, caic, icl, ll)
#    Winner = apply(out, 1, function(x) (1:length(x))[x == 
#        max(x)])
    win = apply(out, 1, which.max)
    rownames(out) = c("AIC", "BIC", "CAIC", "ICL", "Loglik")
    colnames(out) = sort(comps)
    Winner = as.numeric(colnames(out)[win])
    cbind(out, Winner)
}


