rigidTransform <- function(theta, p0, N, cen) {

    if(missing(N)) N <- dim(p0)[1]

    # cdim <- dim(cen)

    # if(is.null(cdim) && length(cen) == 2) cen <- matrix(cen, N, 2, byrow = TRUE)
    # else if(!is.null(cdim)) {

 #      if((cdim[1] == 1 && cdim[2] == 2) || (cdim[1] == 2 && cdim[2] == 1)) cen <- matrix(cen, N, 2, byrow = TRUE)
 #      else if(cdim[1] == 2 && cdim[2] == N) cen <- t(cen)
 #      else if(!(cdim[1] == N && cdim[2] == 2)) stop("rigidTransform: invalid cen argument.")

  #   } else stop("rigidTransform: invalid cen argument.")

#    if(length(theta) %in% c(2, 3)) {
#
#	p0cen <- p0cen - matrix(theta[ 1:2 ], N, 2, byrow = TRUE)
#
#    } else if(length(theta) == 1)

    p0cen <- p0 - cen

#    else stop("rigidTransform: invalid theta argument.")

    # rotation
    r <- theta[ 3 ]
    p0cen <- cbind(cos(r) * p0cen[, 1] - sin(r) * p0cen[, 2], sin(r) * p0cen[, 1] + cos(r) * p0cen[, 2])

    # Phi <- rbind(c(cos(theta[ 3 ]), -sin(theta[ 3 ])),
    #                c(sin(theta[ 3 ]), cos(theta[ 3 ])))
    # p0cen <- t(Phi %*% t(p0cen))

    # translate.
    res <- p0cen + matrix(theta[1:2], N, 2, byrow = TRUE) + cen

    return(res)

} # end of 'rigidTransform' function.

rigider <- function(x1, x0, p0, init = c(0, 0, 0), type = c("regular", "fast"), translate = TRUE, rotate = FALSE, loss, loss.args = NULL,
    interp = "bicubic", method = "BFGS", stages = TRUE, verbose = FALSE, ...) {

    # Below is a little awkward because stages is TRUE by default,
    # which is good if rotate is chosen to be TRUE, but bad otherwise.

    if(stages && !(translate & rotate)) stages <- FALSE

    if(!translate && !rotate) stop("rigider: one of translate or rotate must be true.")

    if(verbose) begin.tiid <- Sys.time()

    theCall <- match.call()

    out <- list()
    out$call <- theCall

    if(missing(x1)) stop("rigider: x1 argument is missing.")
    if(missing(x0)) stop("rigider: x0 argument is missing.")
    if(missing(p0)) stop("rigider: p0 argument is missing.")

    if(missing(loss)) loss <- "QlossRigid"

    type <- tolower(type)
    type <- match.arg(type)

    bigN <- dim(p0)[ 1 ]

    # field.center <- matrix(colMeans(p0, na.rm = TRUE), bigN, 2, byrow = TRUE)
    field.center <- matrix(colMeans(p0[x1 > 0,], na.rm = TRUE), bigN, 2, byrow = TRUE)

    if(type == "regular") {
    
        xdim <- dim(x1)
        if(any(xdim != dim(x0))) stop("rigider: x1 and x0 must have the same dimension.")
    
        outpar <- numeric(0)
    
        if(missing(init)) {
    
            if(verbose) cat("Initial values not passed.  Determining good starting values now.\n")
    
            hold1 <- imomenter(x1, loc = p0)
            hold0 <- imomenter(x0, loc = p0)
    
            tr <- hold1$centroid - hold0$centroid
    
            if(rotate) init <- c(tr, hold1$orientation.angle - hold0$orientation.angle)
            else init <- c(tr, 0)
    
            if(verbose) {
    
                cat( "initial values:\n" )
                print( init )
    
            } # end of if 'verbose' stmt.
    
        } # end of if missing init stmt.
    
        if(!stages && (!translate || !rotate)) stages <- TRUE
    
        if(stages) {
    
            if(translate) {
       
                ofun <- function(theta, p0, x1, x0, loss, loss.args, interp, N, cen, xdim, tr, ...) {
    
    		if(missing(tr)) tr <- NULL
       
                    if(any(c(theta[ 1 ] > xdim[ 1 ] / 2 + 1, theta[ 2 ] > xdim[ 2 ] / 2 + 1))) return(1e16)
       
                    # rigidly transformed coordinates.
                    p1 <- rigidTransform(theta = c(theta, 0), p0 = p0, N = N, cen = cen)
       
                    # rigidly transformed image of x1.
                    y1 <- Fint2d(X = x1, Ws = p1, s = p0, method = interp)
    
                    # loss
                    res <- do.call(loss, c(list(y1 = y1, x0 = x0, p1 = p1, p0 = p0), loss.args))
    
                    return(res)
    
                } # end of internal 'ofun' function.
    
    
                if(verbose) cat("Optimizing translation.\n")
    
                res <- optim(init[ 1:2 ], ofun, p0 = p0, x1 = x1, x0 = x0, loss = loss, loss.args = loss.args,
                    interp = interp, N = bigN, cen = field.center, xdim = xdim, method = method, ...)
    
                if(verbose) cat("Optimal translation found to be: ", res$par, "\nwhere loss value is: ", res$value, "\n")
    
                res$method <- method
    
                p1 <- rigidTransform(theta = c(res$par, init[ 3 ]), p0 = p0, N = bigN, cen = field.center)
                y1 <- Fint2d(X = x1, Ws = p1, s = p0, method = interp)
    
                res$p1 <- p1
                out$x1.translated <- y1
    
                outpar <- res$par
                outval <- res$value
    
            } # end of if 'translate' stmt.
    
    	if(rotate) {
    
                if(translate) init2 <- c(res$par, init[ 3 ])
                else init2 <- init
    
                ofun <- function(theta, p0, x1, x0, loss, loss.args, interp, N, cen, xdim, tr, ...) {
    
                   if((theta > pi / 2) || (theta < - pi / 2)) return(1e16)
    
                   # rigidly transformed coordinates.
                   p1 <- rigidTransform(theta = c(tr, theta), p0 = p0, N = N, cen = cen)
    
                   # rigidly transformed image of x1.
                   y1 <- Fint2d(X = x1, Ws = p1, s = p0, method = interp)
    
                   # loss
                   res <- do.call(loss, c(list(y1 = y1, x0 = x0, p1 = p1, p0 = p0), loss.args))
    
                   return(res)
    
                } # end of internal 'ofun' function.
    
                if(!is.element(method, c("BFGS", "Brent"))) {
    
                    warning("rigider: optimization method is set to something other than Brent or BFGS.  Changing to BFGS for rotation-only.")
                    method2 <- "BFGS"
    
                } else method2 <- method
    
                if(verbose) cat("Optimizing rotation.\n")
    
                res2 <- optim(init2[ 3 ], ofun, p0 = p0, x1 = x1, x0 = x0, loss = loss, loss.args = loss.args,
                            interp = interp, N = bigN, cen = field.center, xdim = xdim, tr = init2[ 1:2 ], method = method2, ...)
    
                if(verbose) cat("Optimal rotation found to be: ", res2$par, "\nwhere loss value is: ", res2$value, "\n")
    
                p1 <- rigidTransform(theta = c(init2[ 1:2 ], res2$par), p0 = p0, N = bigN, cen = field.center)
                y1 <- Fint2d(X = x1, Ws = p1, s = p0, method = interp)
    
                outpar <- c(outpar, res2$par)
                outval <- res2$value
    
            } # end of if 'rotate' stmt.
    
        } else {
    
            ofun <- function(theta, p0, x1, x0, loss, loss.args, interp, N, cen, xdim, tr, ...) {
    
    	    if(missing(tr)) tr <- NULL
    
                if(any(c(theta[ 1 ] > xdim[ 1 ] / 2 + 1, theta[ 2 ] > xdim[ 2 ] / 2 + 1))) return(1e16)
                if((theta > pi / 2) || (theta < - pi / 2)) return(1e16)
    
                # rigidly transformed coordinates.
                   p1 <- rigidTransform(theta = theta, p0 = p0, N = N, cen = cen)
    
                # rigidly transformed image of x1.
                y1 <- Fint2d(X = x1, Ws = p1, s = p0, method = interp)
    
                # loss
                res <- do.call(loss, c(list(y1 = y1, x0 = x0, p1 = p1, p0 = p0), loss.args))
    
                return(res)
    
            } # end of internal 'ofun' function.
    
            if(verbose) cat("Optimizing rigid transformation.\n")
    
            res <- optim(init, ofun, p0 = p0, x1 = x1, x0 = x0, loss = loss, loss.args = loss.args,
                    interp = interp, N = bigN, cen = field.center, xdim = xdim, method = method, ...)
    
            if(verbose) cat("Optimal transformation found to be: ", res$par, "\nwhere loss value is: ", res$value, "\n")
    
            p1 <- rigidTransform(theta = res$par, p0 = p0, N = bigN, cen = field.center)
    
            outpar <- res$par
            outval <- res$val
    
        } # end of if else 'stages' stmts.
    
        if(translate && !rotate) names(outpar) <- c("x", "y")
        else if(!translate && rotate) names(outpar) <- "angle"
        else names(outpar) <- c("x", "y", "angle")
    
        if(stages) {
    
            if(translate) out$translation.only <- res
            if(rotate) out$rotate <- res2
    
        } else out$optim.object <- res
    
    } else if(type == "fast") {

	hold1 <- imomenter(x1, loc = p0)
        hold0 <- imomenter(x0, loc = p0)

        if(translate) tr <- hold1$centroid - hold0$centroid
	if(rotate) rot <- hold1$orientation.angle - hold0$orientation.angle

        if(translate && rotate) {

	    outpar <- inpar <- c(tr, rot)
	    names(outpar) <- c("x", "y", "theta")

	} else if(!translate && rotate) {

	    outpar <- rot
	    inpar <- c(0, 0, rot)
	    names(outpar) <- "theta"

	} else if(translate && !rotate) {

	    outpar <- tr
	    names(outpar) <- c("x", "y")
	    inpar <- c(tr, 0)

	}

	if(stages) {

	    p1.tr <- rigidTransform(theta = tr, p0 = p0, N = bigN, cen = field.center)
	    y1.tr <- Fint2d(X = x1, Ws = p1.tr, s = p0, method = interp)
	    out$x1.translated <- y1.tr

	}

	p1 <- rigidTransform(theta = inpar, p0 = p0, N = bigN, cen = field.center)

    } # end of if 'type' stmts.

    y1 <- Fint2d(X = x1, Ws = p1, s = p0, method = interp)

    out$type <- type

    if(type == "regular") {

	out$initial <- init
	out$value <- outval
	out$optim.args <- list(...)
	out$loss <- list(name = loss, args = loss.args)

    }

    out$interp.method <- interp

    out$par <- outpar
    out$x0 <- x0
    out$x1 <- x1
    out$p0 <- p0
    out$p1 <- p1
    out$x1.transformed <- y1

    class(out) <- "rigided"

    if(verbose) print(Sys.time() - begin.tiid)

    return(out)

} # end of 'rigider' function.

plot.rigided <- function(x, ...) {

    zl <- range(c(c(x$x0), c(x$x1)), finite = TRUE)

    par(mfrow = c(2, 3))
    image(x$x0, col = c("gray", tim.colors(64)), zlim = zl, main = "x0", ...)
    image(x$x1, col = c("gray", tim.colors(64)), zlim = zl, main = "x1", ...)
    image.plot(x$x1.transformed, col = c("gray", tim.colors(64)), zlim = zl, main = "x1 Rigidly Transformed", ...)
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    image.plot(x$x1 - x$x0, col = c("gray", tim.colors(64)), main = "x1 - x0", ...)
    image.plot(x$x1.transformed - x$x0, col = c("gray", tim.colors(64)), main = "x1 rigidly transformed - x0", ...)

    invisible()

} # end of 'plot.rigided' function.


print.rigided <- function(x, ...) {

    print(x$call)

    cat("\n\n")
    cat("Optimal rigid transformation: ", x$par, "\n")
    cat("Objective function value = ", x$value, "\n")

    cat("\n\n")
    if(!is.null(x$translation.only)) cat("Translation convergence code (see ?optim for details): ", x$translation.only$convergence, "\n")
    if(!is.null(x$rotate)) cat("Rotation convergence code (see ?optim for details): ", x$rotate$convergence, "\n")
    if(!is.null(x$optim.object)) cat("Convergence code (see ?optim for details): ", x$optim.object$convergence, "\n")

    invisible()

} # end of 'print.rigided' function.


summary.rigided <- function(object, ...) {

    print(object)
    cat("\n\n")

    MSE0 <- sum(colSums((object$x1 - object$x0)^2, na.rm = TRUE), na.rm = TRUE)
    MSE1 <- sum(colSums((object$x1.transformed - object$x0)^2, na.rm = TRUE), na.rm = TRUE)

    res <- ((MSE0 - MSE1) / MSE0) * 100

    out <- c(MSE0, MSE1, res)
    names(out) <- c("MSE0", "MSE1", "% error reduction")

    print(out)

    invisible(out)

} # end of 'summary.rigided' function.


