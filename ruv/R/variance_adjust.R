variance_adjust <-
function (fit, ebayes = TRUE, evar = TRUE, rsvar = TRUE, bin = 10, 
    rescaleconst = NULL) 
{
    n = ncol(fit$betahat)
    p = nrow(fit$betahat)
    if (TRUE) {
        varbetahat = p.BH = matrix(0, p, n)
        for (l in 1:p) {
            multiplier = fit$multiplier[l]
            varbetahat[l, ] = fit$sigma2 * multiplier
            p.BH[l, ] = p.adjust(fit$p[l, ], method = "BH")
        }
        fit$p.BH = p.BH
        fit$varbetahat = varbetahat
    }
    if (rsvar) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            multiplier = mean(fit$betahat[l, fit$ctl]^2/fit$sigma2[fit$ctl])
            varbetahat[l, ] = fit$sigma2 * multiplier
            tvals[l, ] = fit$betahat[l, ]/sqrt(varbetahat[l, 
                ])
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.rsvar = pvals
        fit$p.rsvar.BH = p.BH
        fit$varbetahat.rsvar = varbetahat
    }
    if (evar) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            varbetahat[l, ] = get_empirical_variances(fit$sigma2, 
                fit$betahat[l, ], bin = bin, rescaleconst = rescaleconst)
            tvals[l, ] = fit$betahat[l, ]/sqrt(varbetahat[l, 
                ])
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), Inf)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.evar = pvals
        fit$p.evar.BH = p.BH
        fit$varbetahat.evar = varbetahat
    }
    if (ebayes) {
        temp = sigmashrink(fit$sigma2, fit$df)
        fit$sigma2.ebayes = temp$sigma2
        fit$df.ebayes = temp$df
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            varbetahat[l, ] = fit$sigma2.ebayes * fit$multiplier[l]
            tvals[l, ] = fit$betahat[l, ]/sqrt(varbetahat[l, 
                ])
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df.ebayes)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.ebayes = pvals
        fit$p.ebayes.BH = p.BH
        fit$varbetahat.ebayes = varbetahat
    }
    if (rsvar & ebayes) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            multiplier = mean(fit$betahat[l, fit$ctl]^2/fit$sigma2.ebayes[fit$ctl])
            varbetahat[l, ] = fit$sigma2.ebayes * multiplier
            tvals[l, ] = fit$betahat[l, ]/sqrt(fit$sigma2.ebayes * 
                multiplier)
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df.ebayes)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.rsvar.ebayes = pvals
        fit$p.rsvar.ebayes.BH = p.BH
        fit$varbetahat.rsvar.ebayes = varbetahat
    }
    return(fit)
}
