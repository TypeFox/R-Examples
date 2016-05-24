lprq <- function(x, y, h, tau = .5, m = 50)
{
	## A toy routine to do locally polynomial quantile regression
        xx <- seq(min(x),max(x),length=m)
        fv <- xx
        dv <- xx
        for(i in 1:length(xx)) {
                z <- x - xx[i]
                wx <- dnorm(z/h)
                r <- rq(y~z, weights=wx, tau=tau, ci=FALSE)
                fv[i] <- r$coef[1.]
                dv[i] <- r$coef[2.]
        }
        list(xx = xx, fv = fv, dv = dv)
}

