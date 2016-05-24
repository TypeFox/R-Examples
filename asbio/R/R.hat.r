R.hat <-  function (M, burn.in = 0.5) {
        alpha = .05                     # 95% intervals
        m = ncol(M)
        x = M [round(((burn.in * nrow(M)) + 1), 0):nrow(M),]  # second half of simulated sequences  h <- round(((burn.in * n) + 1), 0)
        n = nrow(x)
        xdot = as.vector(apply(x, 2, mean))
        s2 = as.vector(apply(x, 2, var))
        W = mean(s2)
        B = n*var(xdot)
        muhat = mean(xdot)
        varW = var(s2)/m
        varB = B^2 * 2/(m-1)
        covWB = (n/m)*(cov(s2,xdot^2) - 2*muhat*cov(s2,xdot))
        sig2hat = ((n-1)*W + B)/n
        quantiles = quantile (as.vector(x), probs=c(.025,.25,.5,.75,.975))

    if (W > 1.e-8) {            # non-degenerate case
        postvar = sig2hat + B/(m*n)
        varpostvar =
                (((n-1)^2)*varW + (1+1/m)^2*varB + 2*(n-1)*(1+1/m)*covWB)/n^2

        chisqdf <- function(A, varA) 2 * (A^2/varA)

        post.df = chisqdf (postvar, varpostvar)
        post.range = muhat + sqrt(postvar) * qt(1-alpha/2, post.df)*c(-1,0,1)
        varlo.df = chisqdf (W, varW)
        confshrink = sqrt(postvar/W*(post.df+3)/(post.df+1))
        confshrink
    }
    else {      # degenerate case:  all entries in "data matrix" are identical
        1
    }
}