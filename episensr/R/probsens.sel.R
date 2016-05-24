#' Probabilistic sensitivity analysis for selection bias.
#'
#' Probabilistic sensitivity analysis to correct for selection bias.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param reps Number of replications to run.
#' @param or.parms List defining the selection bias odds. The first argument provides the probability distribution function (constant, uniform, triangular, or trapezoidal) and the second its parameters as a vector:
#' \enumerate{
#' \item Constant: constant value,
#' \item Uniform: min, max,
#' \item Triangular: lower limit, upper limit, mode,
#' \item Trapezoidal: min, lower mode, upper mode, max.
#' \item Logit-logistic: location, scale, lower bound shift, upper bound shift,
#' \item Logit-normal: location, scale, lower bound shift, upper bound shift.
#' }
#' @param alpha Significance level.
#' @param dec Number of decimals in the printout.
#' @param print A logical scalar. Should the results be printed?
#'
#' @return A list with elements:
#' \item{obs.data}{The analysed 2 x 2 table from the observed data.}
#' \item{obs.measures}{A table of observed odds ratio with confidence intervals.}
#' \item{adj.measures}{A table of corrected odds ratios.}
#' \item{sim.df}{Data frame of random parameters and computed values.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#'
#' @examples
#' # The data for this example come from:
#' # Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N. et al.
#' # Population-based incidence estimates of uveal melanoma in Germany.
#' # Supplementing cancer registry data by case-control data.
#' # Eur J Cancer Prev 2006;15:165-70.
#' set.seed(123)
#' probsens.sel(matrix(c(136, 107, 297, 165),
#' dimnames = list(c("Melanoma+", "Melanoma-"), c("Mobile+", "Mobile-")), nrow = 2, byrow = TRUE),
#' reps = 20000,
#' or.parms = list("triangular", c(.35, 1.1, .43)))
#' @export
#' @importFrom stats median qnorm quantile runif
probsens.sel <- function(case,
                         exposed,
                         reps = 1000,
                         or.parms = list(dist = c("constant", "uniform", "triangular",
                                             "trapezoidal", "logit-logistic",
                                                  "logit-normal"),
                             parms = NULL),
                         alpha = 0.05,
                         dec = 4,
                         print = TRUE){
    if(reps < 1)
        stop(paste("Invalid argument: reps =", reps))
    
    if(is.null(or.parms))
        stop('Please provide odds ratio for the probability of being selected.')
    if(!is.list(or.parms))
        stop('Odds ratio for the probability of being selected should be a list.')
    else or.parms <- or.parms
    if(or.parms[[1]] == "constant" & length(or.parms[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(or.parms[[1]] == "uniform" & length(or.parms[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(or.parms[[1]] == "uniform" & or.parms[[2]][1] >= or.parms[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(or.parms[[1]] == "triangular" & length(or.parms[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(or.parms[[1]] == "triangular" & ((or.parms[[2]][1] > or.parms[[2]][3]) |
                                        (or.parms[[2]][2] < or.parms[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(or.parms[[1]] == "trapezoidal" & length(or.parms[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(or.parms[[1]] == "trapezoidal" & ((or.parms[[2]][1] > or.parms[[2]][2]) |
                                         (or.parms[[2]][2] > or.parms[[2]][3]) |
                                         (or.parms[[2]][3] > or.parms[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')    
        if(or.parms[[1]] == "logit-logistic" & (length(or.parms[[2]]) < 2 | length(or.parms[[2]]) == 3 | length(or.parms[[2]]) > 4))
        stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(or.parms[[1]] == "logit-logistic" & length(or.parms[[2]]) == 4 &
       ((or.parms[[2]][3] >= or.parms[[2]][4]) | (!all(or.parms[[2]][3:4] >= 0 & or.parms[[2]][3:4] <= 1))))
        stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(or.parms[[1]] == "logit-logistic" & length(or.parms[[2]]) == 2)
        or.parms <- list(or.parms[[1]], c(or.parms[[2]], c(0, 1)))
    if(or.parms[[1]] == "logit-normal" & (length(or.parms[[2]]) < 2 | length(or.parms[[2]]) == 3 | length(or.parms[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(or.parms[[1]] == "logit-normal" & length(or.parms[[2]]) == 4 &
       ((or.parms[[2]][3] >= or.parms[[2]][4]) | (!all(or.parms[[2]][3:4] >= 0 & or.parms[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(or.parms[[1]] == "logit-normal" & length(or.parms[[2]]) == 2)
        or.parms <- list(or.parms[[1]], c(or.parms[[2]], c(0, 1)))

    if(inherits(case, c("table", "matrix")))
        tab <- case
    else tab <- table(case, exposed)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    draws <- matrix(NA, nrow = reps, ncol = 4)
    colnames(draws) <- c("or.sel", "corr.or", "reps", "tot.or")

    or.sel <- c(reps, or.parms[[2]])

    obs.or <- (a/b) / (c/d)
    se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
    uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

    logitlog.dstr <- function(sesp) {
        u <- runif(sesp[[1]])
        w <- sesp[[2]] + sesp[[3]] * (log(u / (1 - u)))
        p <- sesp[[4]] + (sesp[[5]] - sesp[[4]]) * exp(w) / (1 + exp(w))
        return(p)
    }
    logitnorm.dstr <- function(sesp) {
        u <- runif(sesp[[1]])
        w <- sesp[[2]] + sesp[[3]] * qnorm(u)
        p <- sesp[[4]] + (sesp[[5]] - sesp[[4]]) * exp(w) / (1 + exp(w))
        return(p)
    }

    if (or.parms[[1]] == "constant") {
        draws[, 1] <- or.parms[[2]]
    }
    if (or.parms[[1]] == "uniform") {
        draws[, 1] <- do.call(runif, as.list(or.sel))
    }
    if (or.parms[[1]] == "triangular") {
        draws[, 1] <- do.call(triangle::rtriangle, as.list(or.sel))
    }
    if (or.parms[[1]] == "trapezoidal") {
        draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(or.sel))
    }
    if (or.parms[[1]] == "logit-logistic") {
        draws[, 1] <- logitlog.dstr(or.sel)
    }
    if (or.parms[[1]] == "logit-normal") {
        draws[, 1] <- logitnorm.dstr(or.sel)
    }

    draws[, 3] <- runif(reps)

    draws[, 2] <- obs.or / draws[, 1]

    draws[, 4] <- exp(log(draws[, 2]) -
                               qnorm(draws[, 3]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))

    corr.or <- c(median(draws[, 2], na.rm = TRUE),
                 quantile(draws[, 2], probs = .025, na.rm = TRUE),
                 quantile(draws[, 2], probs = .975, na.rm = TRUE))
    tot.or <- c(median(draws[, 4], na.rm = TRUE),
                quantile(draws[, 4], probs = .025, na.rm = TRUE),
                quantile(draws[, 4], probs = .975, na.rm = TRUE))        

    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    if (print) {
        cat("Observed Data:",
            "\n--------------", 
            "nOutcome   :", rownames(tab)[1],
            "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
        print(round(tab, dec))
        cat("\n")
        }
    rmat <- data.frame(obs.or, lci.obs.or, uci.obs.or)
    rownames(rmat) <- "Observed Odds Ratio:"
    colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                       sep = ""), "interval")
    if (print) {
        cat("Observed Measures of Exposure-Outcome Relationship:",
            "\n-----------------------------------------------------------------------------------\n\n")
        print(round(rmat, dec))
        cat("\n")
        }
    rmatc <- rbind(corr.or, tot.or)
    rownames(rmatc) <- c("Odds Ratio -- systematic error:", "Odds Ratio -- systematic and random error")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    if (print) {
        print(round(rmatc, dec))
        cat("\nBias Parameters:",
            "\n----------------\n\n")
        cat("OR selection:", or.parms[[1]], "(", or.parms[[2]], ")",
            "\n")
        }
    invisible(list(obs.data = tab,
                   obs.measures = rmat, 
                   adj.measures = rmatc, 
                   sim.df = as.data.frame(draws[, -3])))
}
