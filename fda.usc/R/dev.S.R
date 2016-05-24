################################################################################
dev.S=function (y, S, obs,family = gaussian(),off,offdf,criteria="GCV",
	W = diag(1, ncol = ncol(S), nrow = nrow(S)), trim = 0, draw = FALSE, ...)
{
    n = ncol(S)
    mu = family$linkinv(S%*%y+off)
    sigma=sqrt(family$variance(mu))
    tab = list("GCV", "AIC", "FPE", "Shibata", "Rice")
    type.i = pmatch(criteria, tab)
    e = (obs - mu)/sigma
    if (trim > 0) {
            e.trunc = quantile(abs(e), probs = (1 - trim), na.rm = TRUE,
                type = 4)
            l <- which(abs(e) <= e.trunc)
    }
    else {  l = 1:n  }
    res = traza(t(e[l]) %*% W[l,l] %*% e[l])
    ndf <- sum(diag(S)[l],na.rm=TRUE)+offdf
    if (is.na(type.i)) {
        if (ndf>0.8*n)  vv = Inf
        else vv = 1/(1 - 2 * ndf/n)
    }
    else {
        vv<-switch(type.i,
                   "1"=if (ndf>0.8*n){vv=Inf} else {vv = (1 - ndf/n)^(-2)},
                   "2"=if (ndf>0.8*n){vv=Inf} else {vv = exp(2 * ndf/n)},
                   "3"=if (ndf>0.8*n){vv=Inf} else {vv = (1 + ndf/n)/(1 - ndf/n)},
                   "4"=if (ndf>0.8*n){vv=Inf} else {vv = (1 + ndf/n)/(1 - ndf/n)},
                   "5"=if (ndf>0.8*n){vv=Inf} else { vv = 1/(1 - 2 * ndf/n)})
         }
    return(res * vv/n)
}

