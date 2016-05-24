IWLSiteration <- function(x, y, inib, iniscale, maxiter, tol, b, c1, c2) {
    # slightly changed subfunction for the Fast-Tau algorithm for linear regression originally published in 
    # Salibian-Barrera, M., Willems, G. und Zamar, R. (2008): The fast-tau estimator for regression.
    # Journal of Computational and Graphical Statistics, 17 Nr. 3, 659-682.

    n <- nrow(x)
    p <- ncol(x)

    res <- y - x %*% inib
    if (iniscale == 0) {
        scale <- median(abs(res))/.6745
    } else {
        scale <- iniscale
    }
    oldbeta <- inib

    betadiff <- 2*tol
    iter <- 0

    fwOpt <- function(x, cc) {
        tmp <- (-1.944 / cc^2 + 1.728 * x^2 / cc^4 - 0.312 * x^4 / cc^6 + 0.016 * x^6 / cc^8) / 3.25
        tmp[abs(x) < 2*cc] <- 1 / (3.25*cc^2)
        tmp[abs(x) > 3*cc] <- 0
        tmp
    }
	
    psixOpt <- function(x, cc) {
        tmp <- x^2 / (3.25*cc^2)
        tmp2 <- (-1.944 * x^2 / cc^2 + 1.728 * x^4 / cc^4 - 0.312 * x^6 / cc^6 + 0.016 * x^8 / cc^8) / 3.25
        tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
        tmp[abs(x) > 3*cc] <- 0
        tmp
    }

    rhoOpt <- function(x, cc) {
        tmp <- x^2 / 2 / (3.25*cc^2)
        tmp2 <- (1.792 - 0.972 * x^2 / cc^2 + 0.432 * x^4 / cc^4 - 0.052 * x^6 / cc^6 + 0.002 * x^8 / cc^8) / 3.25
        tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
        tmp[abs(x) > 3*cc] <- 1
        tmp
    }

    WtellerOpt <- function(x, cc) {
        tmp <- (3.584 - 0.864 * x^4 / cc^4 + 0.208 * x^6 / cc^6 - 0.012 * x^8 / cc^8) / 3.25
        tmp[abs(x) < 2*cc] <- 0
        tmp[abs(x) > 3*cc] <- 2
        tmp
    }

    while ((betadiff > tol) && (iter < maxiter)) {
        scale <- sqrt( scale^2 * mean( rhoOpt(res/scale,c1) ) / b )
        scaledres <- res/scale
        Wn.teller <- sum(WtellerOpt(scaledres,c2))
        Wn.noemer <- sum(psixOpt(scaledres,c1))
        Wn <- Wn.teller / Wn.noemer
        weights <- (Wn * fwOpt(scaledres,c1) + fwOpt(scaledres,c2))
        weights <- apply(cbind(weights, 0), 1, max) # in order not to get an infenitesemal small, negative number
        sqweights <- sqrt(weights)
        xw <- x * as.vector(sqweights)
        yw <- y * sqweights
        newbeta <- qr.coef(qr(xw),yw)
        if (any(!is.finite(newbeta))) {
            newbeta <- inib
            scale <- iniscale
            break
        }
        betadiff <- sqrt(sum((oldbeta - newbeta)^2))
        res <- y - x %*% newbeta
        oldbeta <- newbeta
        iter <- iter + 1
    }
	
    return( list( betarw = newbeta, scalerw = scale ) )

}
