mwtp <- function(output, monetary.variables, nonmonetary.variables = NULL, nreplications = 10000,
                 percentile.points = NULL, confidence.level = 0.95, method = "kr", seed = NULL) 
{
# Name: mwtp
# Title: Calculating the marginal willingness to pay
# Arguments:
#  output                  An object containing the output from clogit() or glm().
#  monetary.variables      A vector containing the names of the monetary variables
#                            in the output used to calculate the MWTPs.
#  nonmonetary.variables   A vector containing the names of the non-monetary variables 
#                            in the output used to calculate the MWTPs.
#  nreplications           An integer value denoting the number of replications
#                            in the simulation method.
#  percentile.points       It is only kept for giving an error message regarding unused argument.
#                            It will be removed. confidence.level is used instead.
#  confidence.level        A value showing the confidence level (coefficient).
#  method                  A character variable describing the method used for calculating MWTPs:
#                            "kr" if Krinsky and Robb's method is used;
#                            "delta" if Delta method is used.
#  seed                    Seed for a random number generator.



# check arguments

    if (is.null(percentile.points) == FALSE) {
        stop(message = "confidence.level should be used instead of percentile.points")
    }
    if (length(confidence.level) > 1) {
        stop(message = "confidence.level should be a value")
    }
    if (any(confidence.level > 1 | confidence.level < 0) == TRUE) {
        stop(message = "confidence.level should be [0, 1]")
    }
    if (length(monetary.variables) > 1) {
        if (is.list(nonmonetary.variables) == FALSE) {
            stop(message = "nonmonetary.variables should be set in a list format")
        }
    }


# set variables, vectors, and matrix

    est.cov  <- vcov(output)         # estimated covariances with variable names
    est.coef <- coefficients(output) # estimated coefficients with variable names

    lower.limit <- (1 - confidence.level) / 2     # lower limit of confidence interval
    upper.limit <- confidence.level + lower.limit # upper limit of confidence interval
    z           <- qnorm(upper.limit)  

    num.monetary.variables    <- length(monetary.variables) # number of monetary variables

    monetary.index <- rep(0, times = num.monetary.variables) 
    for (i in 1:num.monetary.variables) {
        monetary.index[i] <- which(names(est.coef) == monetary.variables[i])
    }

    if (is.null(nonmonetary.variables) == TRUE) {
        nonmonetary.variables <- names(est.coef)[-monetary.index]
    }

    unlist.nonmonetary.variables <- unlist(nonmonetary.variables)
    num.nonmonetary.variables <- length(unlist.nonmonetary.variables)

    nonmonetary.index <- rep(0, times = num.nonmonetary.variables)
    for (i in 1:num.nonmonetary.variables) {
        nonmonetary.index[i] <- which(names(est.coef) == unlist.nonmonetary.variables[i])
    }

    if (is.list(nonmonetary.variables) == TRUE) { # labeled design
        each.times <- sapply(nonmonetary.variables, length)
        vrep       <- Vectorize(rep)
        monetary.index.expanded <- vrep(monetary.index, each.times)
    }

# calculate confidence intervals

    # Krinsky and Robb's method
    if(method == "kr") {

        # replicate coefficients using Krinsky and Robb's method
        if (is.null(seed) == FALSE) {
            set.seed(seed)
        }
        repb <- mvrnorm(nreplications, est.coef, est.cov)

        # calculate mean MWTPs and simulate empirical distribution for each of MWTPs 
        # labeled design
        if (is.list(nonmonetary.variables) == TRUE) {
            mwtps    <- -est.coef[nonmonetary.index] / est.coef[monetary.index.expanded]
            repmwtps <- -repb[, nonmonetary.index] / repb[, monetary.index.expanded]
        }
        # unlabeled design
        else {
            mwtps    <- -est.coef[nonmonetary.index] / est.coef[monetary.index]
            repmwtps <- -repb[, nonmonetary.index] / repb[, monetary.index]
        }

        # calculate confidence intervals
        confidence.intervals <- apply(repmwtps, 2, quantile, 
                                      probs = c(lower.limit, upper.limit))

        # format output
        rtn <- list(mwtp.table = t(rbind(MWTP = mwtps, confidence.intervals)),
                    method     = method,
                    mwtps      = repmwtps,
                    repb       = repb)

    }

    # Delta method
    else {

        # define function for calculating variance of MWTP
        variance.mwtps <- function(coef.monetary,
                                   var.monetary,
                                   coef.nonmonetary,
                                   var.nonmonetary,
                                   cov.nonm.mon) {
            return((var.nonmonetary / coef.monetary^2) +
                   ((coef.nonmonetary^2 * var.monetary) / coef.monetary^4) -
                   ((2 * coef.nonmonetary * cov.nonm.mon) / coef.monetary^3))
        }

        # calculate MWTPs and their variances
        # labeled design
        if (is.list(nonmonetary.variables) == TRUE) {
            mwtps <- -est.coef[nonmonetary.index] / est.coef[monetary.index.expanded] # MWTP
            # extract covariances needed for calculating variances of MWTPs from estimated covariances
            cov   <- rep(0, num.nonmonetary.variables)
            for (i in 1:num.nonmonetary.variables) {
                cov[i] <- est.cov[nonmonetary.index[i], monetary.index.expanded[i]]
            }
            var.mwtps <- variance.mwtps(coef.monetary    = est.coef[monetary.index.expanded],
                                        var.monetary     = diag(est.cov)[monetary.index.expanded],
                                        coef.nonmonetary = est.coef[nonmonetary.index],
                                        var.nonmonetary  = diag(est.cov)[nonmonetary.index],
                                        cov.nonm.mon     = cov)
        }
        # unlabeled design
        else {
            mwtps     <- -est.coef[nonmonetary.index] / est.coef[monetary.index]
            cov       <- est.cov[nonmonetary.index, monetary.index]
            var.mwtps <- variance.mwtps(coef.monetary    = est.coef[monetary.index],
                                        var.monetary     = diag(est.cov)[monetary.index],
                                        coef.nonmonetary = est.coef[nonmonetary.index],
                                        var.nonmonetary  = diag(est.cov)[nonmonetary.index],
                                        cov.nonm.mon     = cov)
        }

        # calculate confidence intervals
        upper.mwtps <- mwtps + z * sqrt(var.mwtps)
        lower.mwtps <- mwtps - z * sqrt(var.mwtps)

        # format output
        mwtp.table <- cbind(mwtps, lower.mwtps, upper.mwtps)
        colnames(mwtp.table) <- c("MWTP",                                   # name of 1st column
                                  paste(lower.limit * 100, "%", sep = ""),  # name of 2nd column
                                  paste(upper.limit * 100, "%", sep = ""))  # name of 3rd column
        rtn <- list(mwtp.table = mwtp.table,
                    method     = method)
    }

# return output

    class(rtn) <- "mwtp"
  
    return(rtn)
}

