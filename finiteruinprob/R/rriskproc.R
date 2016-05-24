rriskproc <- function(m      = 1001L,
                      window = c(0.0, 1.0),
                      num    = 1L,
                      sigma  = 1.0,
                      freq   = 1.0,
                      drift  = 0.0,
                      jumpdist, ...) {

    stopifnot(length(window) == 2L,
              all(is.finite(window)))

    m          <- as.integer(m)
    num        <- as.integer(num)
    sigma      <- as.double(sigma)
    freq       <- as.double(freq)
    drift      <- as.double(drift)
    timerange  <- range(as.double(window), na.rm = TRUE)
    .frequency <- (m - 1.0) / diff(timerange)
    .time      <- seq.int(timerange[1L], timerange[2L], by = 1.0 / .frequency)


    ## BROWNIAN MOTION
    ## ===============
    ##
    ## Since this method generates two (more generally: an even number) of
    ## independent Brownian motions simultaneously, calculate how many pairs
    ## are needed.
    bmnum <- num %/% 2L + 1L

    ## Simulate multiple Brownian motions with volatility sigma using the
    ## circulant embedding method as provided by
    ##     Dietrich and Newsam (1997).  Fast and exact simulation of
    ##     stationary Gaussian processes through circulant embedding of the
    ##     covariance matrix.  SIAM Journal of Scientific Computing, 18(4),
    ##     pp. 1088--1107.
    ##
    # k1       <- exp(-.time)
    # eigenval <- Re(fft(c(k1, k1[(m - 1L):2L])))
    # eigenval <- Re(fft(exp(-.time)[c(1L:m, (m - 1L):2L)]))
 
    m2 <- 2L * m - 2L

    noise <- matrix(complex(real      = rnorm(bmnum * m2),
                            imaginary = rnorm(bmnum * m2)),
                    ncol = bmnum,
                    nrow = m2)

    # mybm <- mvfft(sweep(noise, 1L, sqrt(eigenval / m2), '*'))
    mybm <- mvfft(sweep(x      = noise,
                        MARGIN = 1L,
                        STATS  = sqrt(Re(fft(exp(-.time)[c(1L:m, (m - 1L):2L)])) / m2),
                        FUN    = '*'))

    ## Subtract the first value of each process (corresponding to time 0) such
    ## that all of the processes start at the origin.
    mybm <- sweep(x      = mybm[1L:m, , drop = FALSE],
                  MARGIN = 2L,
                  STATS  =  mybm[1L, , drop = FALSE],
                  FUN    = '-')


    ## AGGREGATE LOSS PROCESS
    ## ======================
    ##
    ## Use the expected number of claims in the observed time window as
    ## initial guess for the number of claims to be generated (and add 20
    ## percent to err on the safe side).  Then generate that many random claim
    ## arrival times.
    nclaims    <- as.integer(freq * diff(timerange) * 1.2)
    myarrivals <- apply(matrix(data = rexp(nclaims * num, freq),
                               nrow = nclaims,
                               ncol = num), 2L, cumsum)

    # TODO: Empty list elements of myarrivals are currently dropped.
    #       This is undesired and leads to an error.

    myarrivals <- lapply(
        ## Make it a list
        unname(as.list(as.data.frame(myarrivals))),
        function(.arrivals) {
            ## If the number of generated claims was not high enough, add half
            ## of the same amount again.  Then repeat this until the latest
            ## claim occurs outside the observation window and prune the claim
            ## arrival times to the observation window.
            if (length(.arrivals) > 0L) {
                while (max(.arrivals) < timerange[2L]) {
                    .arrivals <- c(.arrivals,
                                   tail(.arrivals, 1L)
                                   + cumsum(rexp(n    = as.integer(ceiling(0.5 * nclaims)),
                                                 rate = freq)))
                }
            }
            .ind <- which(timerange[1L] <= .arrivals & .arrivals <= timerange[2L])
            ## Round the claim arrival times such that they lie on the same
            ## grid as the Brownian motion.
            return(round(.frequency * .arrivals[.ind]) / .frequency)
        }
    )

    ## Generate the cumulated claim amounts corresponding to those claim
    ## arrival times.
    myclaims <- lapply(vapply(myarrivals, length, integer(1L)),
                       function(n) c(0.0, cumsum(jumpdist(n, ...))))

    ## Return the value of the aggregate claim amount process at the grid
    ## points by using a step function with steps at the gridified claim
    ## arrival times.
    claimproc <- vapply(seq_len(num),
                        function(.ind) {
                            approx(x      = c(0.0, myarrivals[[.ind]]),
                                   y      = myclaims[[.ind]],
                                   xout   = .time,
                                   method = 'constant',
                                   rule   = 2L,
                                   f      = 0.0,
                                   ties   = max)$y
                        },
                        double(m))

    ## Put everything together
    .data <- unname(sigma * cbind(Re(mybm), Im(mybm))[, seq_len(num)]
                    + drift * .time - claimproc)

    return(ts(data      = .data,
              start     = timerange[1L],
              end       = timerange[2L],
              frequency = .frequency))
}
