re.s <- function(x,y,initial.beta,initial.scale,kk,conv,b,cc) {
    # slightly modified subfunction for the Fast-S algorithm for linear regression originally published in 
    # Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
    # Journal of Computational and Graphical Statistics, 15 (2), 414-427

    # does "kk" IRWLS refining steps from "initial.beta"
    #
    # if "initial.scale" is present, it's used, o/w the MAD is used
    # kk = number of refining steps
    # conv = 0 means "do kk steps and don't check for convergence"
    # conv = 1 means "stop when convergence is detected, or the
    #                 maximum number of iterations is achieved"
    # b and cc = tuning constants of the equation
    #

    n <- dim(x)[1]
    p <- dim(x)[2]
    res <- y - x %*% initial.beta
    if( missing( initial.scale ) ) initial.scale <- scale <- median(abs(res))/.6745 else scale <- initial.scale

    if( conv == 1) kk <- 50
    #
    # if conv == 1 then set the max no. of iterations to 50
    # magic number alert!!!

    beta <- initial.beta

    f.w <- function(u, cc) {
        # weight function = psi(u)/u
        tmp <- (1 - (u/cc)^2)^2
        tmp[ abs(u/cc) > 1 ] <- 0
        return(tmp)
    }
    our.solve <- function(a,b) {
        a <- qr(a)
        da <- dim(a$qr)
        if(a$rank < (p <- da[2])) return(NA) else qr.coef(a, b)
    }
    rho <- function(u, cc=1.56) {
        w <- abs(u)<=cc
        v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
        v <- v*6/cc^2
        return(v)
    }

    for(i in 1:kk) {
        # do one step of the iterations to solve for the scale
        scale.super.old <- scale
        scale <- sqrt( scale^2 * mean( rho( res / scale, cc ) ) / b     )
        # now do one step of IRWLS with the "improved scale"
        weights <- f.w( res/scale, cc )
        W <- matrix(weights, n, p)
        xw <- x * sqrt(W)
        yw <- y *   sqrt(weights)
        beta.1 <- our.solve( t(xw) %*% xw ,t(xw) %*% yw )
        if(any(is.na(beta.1))) {
            beta.1 <- initial.beta
            scale <- initial.scale
            break
        }
        if( (conv==1) ) {
            # check for convergence
            if( norm( beta - beta.1, "2" ) / norm(beta, "2") < 1e-20 ) break
            # magic number alert!!!
        }
        res <- y - x %*% beta.1
        beta <- beta.1
    }
	
    res <- y - x %*% beta
    # get the residuals from the last beta
    return(list(beta.rw = beta.1, scale.rw = scale))

}
