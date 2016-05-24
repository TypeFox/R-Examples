##
##  n e l m i n . R  Nelder-Mead function minimization
##


nelmin <- function(fn, x0, tol = 1e-10, ..., 
            maxfeval = 1e4, step = rep(1.0, length(x0))) {
    stopifnot(is.numeric(x0), is.numeric(step))
    n <- length(x0)
    if (length(step) != n)
        stop("Argument 'step' must be of the same length as 'x0'.")


    fun <- match.fun(fn)
    fn  <- function(x) fun(x, ...)

                        # Inputs:
    start  <- x0        # starting point
    reqmin <- tol       # the terminating limit for the variance of function values
    # step <- step      # size and shape of the initial simplex
    kcount <- maxfeval  # maximum number of function evaluations.
    konvge <- kcount/100# convergence check is carried out every KONVGE iterations, >= 1

                        # Outputs:
    xmin   <- NA        # estimated minimum of the function
    ynewlo <- NA        # minimum value of the function
    icount <- 0         # number of function evaluations
    numres <- 0         # number of restarts, must be > 1
    ifault <- 0         # error indicator, 0, 1, 2

    # Constants for Nelder-Mead
    ccoeff <- 0.5
    ecoeff <- 2.0
    eps    <- 0.001
    rcoeff <- 1.0
    
    jcount <- konvge
    dn <- n
    nn <- n + 1
    dnn <- nn
    del <- 1.0
    rq <- reqmin * dn

    pbar   <- numeric(n)  # centroid
    pstar  <- numeric(n)
    p2star <- numeric(n)

    y <- numeric(n+1)
    p <- matrix(0, n, n+1)

    while ( TRUE ) {                    # outer while loop
        
        # Initial or restarted loop
        p[, nn] <- start
        y[nn] <- fn ( start )
        icount <- icount + 1

        for (j in 1:n) {
          x <- start[j]
          start[j] <- start[j] + step[j] * del
          p[, j] <- start
          y[j] <- fn ( start )
          icount <- icount + 1
          start[j] <- x
        } # simplex construction is complete

        # Find highest and lowest y values.
        ilo <- which.min(y)
        ylo <- y[ilo]

        while ( icount < kcount ) {                # inner while loop

            # indicate the vertex of the simplex to be replaced
            ihi <- which.max(y)
            ynewlo <- y[ihi]

            # Calculate the centroid of the simplex vertices
            #   excepting the vertex with Y value YNEWLO
            pbar <- rowSums(p[, -ihi]) / dn

            # Reflection through the centroid
            pstar <- pbar + rcoeff * ( pbar - p[,ihi] )
            ystar <- fn ( pstar )
            icount <- icount + 1

            # Successful reflection, so extension
            if ( ystar < ylo ) {
                p2star = pbar + ecoeff * ( pstar - pbar )
                y2star <- fn ( p2star )
                icount <- icount + 1

                # Check extension.
                if ( ystar < y2star ) {
                    p[, ihi] <- pstar
                    y[ihi] <- ystar

                # Retain extension or contraction.
                } else {
                    p[, ihi] <- p2star
                    y[ihi] <- y2star
                }

            #  No extension.
            } else {
                l <- sum(ystar < y)
                if ( l > 1 ) {
                    p[, ihi] <- pstar
                    y[ihi] <- ystar

                # Contraction on the Y(IHI) side of the centroid.
                } else if ( l == 0 ) {
                    p2star <- pbar + ccoeff * ( p[, ihi] - pbar )
                    y2star <- fn ( p2star )
                    icount <- icount + 1

                    #  Contract the whole simplex.
                    if ( y[ihi] < y2star ) {
                        for (j in 1:nn) {
                            p[, j] <- (p[, j] + p[, ilo]) / 2
                            xmin <- p[, j]
                            y[j] <- fn ( xmin )
                            icount <- icount + 1
                        }

                        ilo <- which.min(y)
                        ylo <- y[ilo]
                        next

                    # Retain contraction
                    } else {
                        p[, ihi] <- p2star
                        y[ihi] <- y2star
                    }

                #  Contraction on the reflection side of the centroid
                } else if ( l == 1 ) {
                    p2star <- pbar + ccoeff * ( pstar - pbar )
                    y2star <- fn ( p2star )
                    icount <- icount + 1

                    #  Retain reflection?
                    if ( y2star <= ystar ) {
                        p[, ihi] <- p2star
                        y[ihi] <- y2star

                    } else {
                        p[, ihi] <- pstar
                        y[ihi] <- ystar
                    }
                }
            }

            #  Check if YLO improved.
            if ( y[ihi] < ylo ) {
                ylo <- y[ihi]
                ilo <- ihi
            }
            jcount <- jcount - 1
            
            if ( 0 < jcount )
                next

            #  Check to see if minimum reached.
            if ( icount <= kcount ) {
              jcount <- konvge

              x <- sum(y) / dnn
              z <- sum((y - x)^2)
              if ( z <= rq )
                break
            }
        }                               # end inner while loop

        # Factorial tests to check that YNEWLO is a local minimum
        xmin <- p[, ilo]
        ynewlo <- y[ilo]

        if ( kcount < icount ) {
            ifault <- 2
            break
        }

        ifault <- 0

        # Check in all directions with step length
        for (i in 1:n) {
            del <- step[i] * eps
            xmin[i] <- xmin[i] + del
            z <- fn ( xmin )
            icount <- icount + 1
            if ( z < ynewlo ) {
                ifault <- 2
                break
            }
            xmin[i] <- xmin[i] - del - del
            z <- fn ( xmin )
            icount <- icount + 1
            if ( z < ynewlo ) {
                ifault <- 2
                break
            }
            xmin[i] <- xmin[i] + del
        }

        if ( ifault == 0 )
            break

        #  Restart the procedure.
        start <- xmin

        del <- eps
        numres <- numres + 1
    }                                   # end outer while loop

    return(list(xmin = xmin, fmin = ynewlo, fcount = icount, restarts = numres))
} # end of function


nelminb <- function (fn, x0, lower, upper, tol = 1e-10, ..., 
                     maxfeval = 10000, step = rep(1, length(x0))) {
	Trf <- transfinite(lower, upper, length(x0))
	h <- Trf$h; hinv <- Trf$hinv

	f <- function(x) fn(hinv(x), ...)  # f must be defined on all of R^n
    S <- nelmin(f, h(x0), tol = tol, maxfeval = maxfeval, step = step)
    S$xmin <- hinv(S$xmin)

    return(S)
}
