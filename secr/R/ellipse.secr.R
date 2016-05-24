############################################################################################
## package 'secr'
## ellipse.secr.R
## confidence ellipse for two named beta parameters
## bivariate normal confidence ellipse for centre of 2-D distribution
## last changed 2010 02 15
############################################################################################

ellipse.secr <- function (object, par=c('g0', 'sigma'), alpha = 0.05, npts = 100,
    plot = TRUE, linkscale = TRUE, add = FALSE, col = palette(), ...) {

    if (inherits(object, 'list')) {
        temp <- list()
        nsecr <- length(object)
        colr <- rep(col, nsecr)
        for (i in 1:nsecr) {
            temp[[i]] <- ellipse.secr (object[[i]], par = par, alpha = alpha, npts = npts,
                plot = plot, linkscale = linkscale, add = add, col = colr[i], ...)
        }
        invisible(temp)
    }
    else {
        if ((length(par) != 2) | (any(is.na(match(par, object$betanames)))))
            stop ("require two named beta parameters")
        if ((!linkscale) & (any(is.na(match(par, object$realnames)))))
            stop ("'ellipse.secr' is currently limited to null models")
        vcv <- vcov(object)[par,par]
        centre <- coef(object)[par,'beta']
        scale  <- sqrt(diag(vcv))
        r <- vcv[par[1],par[2]] / prod(scale)
        d <- acos(r)
        X2 <- qchisq(1-alpha, 2)
        t <- sqrt(X2)
        pairwise <-  solve(vcv)
        a <- 2 * pi * (0:npts) / npts
        x <- centre[1] + t * scale[1] * cos(a + d/2)
        y <- centre[2] + t * scale[2] * cos(a - d/2)
        if (!linkscale) {
            x <- untransform (x, object$link[[par[1]]])
            y <- untransform (y, object$link[[par[2]]])
            centre[1] <- untransform (centre[1], object$link[[par[1]]])
            centre[2] <- untransform (centre[2], object$link[[par[2]]])
        }
        if (plot) {
            if (add)
                lines(x, y, col = col, ...)
            else
                plot(x,y,type='l', xlab=par[1], ylab=par[2], col = col, ...)
                points (centre[1], centre[2], pch=3, col = col)
        }
        invisible(list(x=x,y=y))
    }
}

## data(secrdemo)
## ellipse.secr(secrdemo.0)
## ellipse.secr(secrdemo.0, link = FALSE, xlim = c(0.1,0.4),
##     ylim = c(20,40), alpha=0.01)
#############################################################################################

ellipse.bvn <- function (xy, alpha = 0.05, npts = 100, centroid = TRUE,
                                add = FALSE, ...) {
    ## xy is 2-column matrix of coordinates for one animal
    centre <- matrix(apply(xy, 2, mean), ncol = 2)
    vcv <- var(xy)
    if (centroid) vcv <- vcv/nrow(xy)
    detS <- det(vcv)^0.5  ## sqrt(generalised variance)
    scale  <- sqrt(diag(vcv))
    r <- vcv[1,2] / prod(scale)
    d <- acos(r)
    X2 <- qchisq(1-alpha, 2)
    t <- sqrt(X2)
    pairwise <-  solve(vcv)
    a <- 2 * pi * (0:npts) / npts
    x <- centre[1] + t * scale[1] * cos(a + d/2)
    y <- centre[2] + t * scale[2] * cos(a - d/2)
    if (!add)
        plot(x,y, type='n', xlab='', ylab='')
    polygon(x, y, ...)
    out <- list(x = x, y = y, centre = centre)
    invisible(out)
}
#############################################################################################
