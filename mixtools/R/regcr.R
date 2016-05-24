# Produce credible region for regression lines based on
# sample from posterior distribution of beta parameters.
# It is assumed that beta, upon entry, is an nx2 matrix,
# where the first column gives the intercepts and the
# second column gives the slopes.
# alpha is the proportion of beta to remove from the posterior.
# Thus, 1-alpha is the level of the credible region.
#
# If nonparametric=TRUE, then the region is based on the convex
# hull of the remaining beta after trimming, which is accomplished
# using a data depth technique.
#
# If nonparametric=FALSE, then the region is based on the
# asymptotic normal approximation.

regcr=function (beta, x, em.beta=NULL, em.sigma=NULL, alpha = 0.05, nonparametric = FALSE, plot = FALSE, 
    xyaxes = TRUE, ...) 
{
    if (nonparametric) {
	  beta.old=beta
	  if(is.null(em.sigma) && is.null(em.beta)){
	  beta=t(matsqrt(solve(cov(beta)))%*%(t(beta)-apply(beta,2,mean)))
	  } else beta=1/sqrt(length(x))*t(em.sigma^(-1)*matsqrt(t(cbind(1,x))%*%cbind(1,x))%*%(t(beta.old)-em.beta))
        d = depth(beta, beta)
        beta = beta.old[order(d), ]
        d = d[order(d)]
        n = length(d)
        trimbeta = beta[-(1:round(n * alpha)), ]
        h = unique(trimbeta[chull(trimbeta), ])
        nh = nrow(h)
        m = which.max(h[, 2])
        h = rbind(h[m:nh, ], h[((1:nh) < m), ], h[m, ])
        bound = NULL
        for (i in 1:nh) {
            bound = rbind(bound, h[i, ])
            bound = rbind(bound, cbind(seq(h[i, 1], h[i + 1, 
                1], len = 50), seq(h[i, 2], h[i + 1, 2], len = 50)))
        }
	beta=trimbeta
    }
    else {
        xbar = apply(beta, 2, mean)
        n = nrow(beta)
        cbeta = t(t(beta) - xbar)
        S = t(cbeta) %*% cbeta/(n - 1)
        eS = eigen(S)
        B = eS$vec %*% diag(sqrt(eS$val))
        theta = seq(0, 2 * pi, len = 250)
        v = cbind(cos(theta), sin(theta)) * sqrt(qchisq(1 - alpha, 
            2))
        h = t(B %*% t(v) + xbar)
        nh = nrow(h)
        m = which.max(h[, 2])
        h = rbind(h[m:nh, ], h[((1:nh) < m), ], h[m, ])
        bound = h
    }
    z <- length(x) * 5
    u <- seq(min(x), max(x), length = z)
    lower <- c()
    upper <- c()
    v <- c()
    for (j in 1:z) {
        for (i in 1:nrow(beta)) {
            v[i] <- as.matrix(beta[, 1][i] + beta[, 2][i] * u[j])
        }
        uv <- cbind(u[j], v)
        lower <- rbind(lower, uv[order(v), ][1, ])
        upper <- rbind(upper, uv[order(v), ][nrow(beta), ])
    }
    if (plot) {
        if (xyaxes) {
            lines(upper, ...)
            lines(lower, ...)
        }
        else {
            lines(bound, ...)
        }
    }
    invisible(list(boundary = bound, upper = upper, lower = lower))
}



