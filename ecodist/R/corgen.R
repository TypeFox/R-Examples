corgen <- function (len, x, r, population = FALSE, epsilon = 0) 
{
# remember the sign of the correlation for later
    rsign <- sign(r)
    if (rsign == 0) 
        rsign <- 1
    r <- abs(r)

# need to either specify x-data or the desired vector length 
    if (missing(x)) {
        if (missing(len)) {
            stop("Must specify x or len.\n")
        }
        else {
            # if x wasn't given, sample it from a normal distribution
            x.rand <- TRUE
            x <- scale(rnorm(len))
            x.orig <- x
        }
    }
    else {
    # can either draw from a population or simulate an exact 
    # correlation (within epsilon), but not both
        if (population == TRUE) {
            # if x is given, population is set to FALSE
            warning("If x is specified population is ignored.\n")
            population <- FALSE
        }
        x.rand <- FALSE
        len <- length(x)
        x.orig <- x
        x <- scale(x)
    }
    if (epsilon != 0) {
        if (population == TRUE) {
            warning("If epsilon is specified population is ignored.\n")
            population <- FALSE
        }
    }

## Here's where the real work starts
    # First, draw y from a normal distribution
    y <- scale(rnorm(len))
    if (!population & x.rand) {
        # if exact correlations are desired, use princomp to
        # create uncorrelated data
        # UNLESS x is given - princomp would trash the
        # specified x values
        xy <- princomp(cbind(x, y))$scores
        x <- xy[, 1]
        x.orig <- x
        y <- xy[, 2]
    }

    # create a new y based on x and the desired correlation
    a <- r/sqrt(1 - r^2)
    y <- x * a + y

    if (epsilon > 0) {
        itcounter <- 0
        # check to see if cor(x, y) meets the given tolerances
        # this is kludgy, but works
        while (abs(cor(x, y) - r) > epsilon) {
            # sometimes takes too long or doesn't converge
            # with given starting point
            if(itcounter > 1000) {
                if(x.rand) {
                    x <- scale(rnorm(len))
                    x.orig <- x
                    itcounter <- 0
                }
                else {
                    # hopefully never get here
                    stop("Doesn't converge.\n")
                }
            }
            itcounter <- itcounter + 1       
            # if not within epsilon, generate a new y
            y <- scale(rnorm(len))
            if (!population & x.rand) {
                xy <- princomp(cbind(x, y))$scores
                x <- xy[, 1]
                x.orig <- x
                y <- xy[, 2]
            }
            a <- r/sqrt(1 - r^2)
            y <- x * a + y
        }
    }

    # return x and y with the correct sign of r restored
    list(x = x.orig, y = y * rsign)
}

