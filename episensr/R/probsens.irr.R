#' Probabilistic sensitivity analysis for exposure misclassification of person-time data and random error.
#'
#' Probabilistic sensitivity analysis to correct for exposure misclassification when person-time data has been collected.
#'
#' @param counts A table or matrix where first row contains disease counts and second row contains person-time at risk, and first and second columns are exposed and unexposed observations, as:
#' \tabular{lll}{
#' \tab Exposed \tab Unexposed \cr
#' Cases \tab a \tab b \cr
#' Person-time \tab N1 \tab N0
#' }
#' @param pt A numeric vector of person-time at risk. If provided, \code{counts} must be a numeric vector of disease counts.
#' @param reps Number of replications to run.
#' @param seca.parms List defining the sensitivity of exposure classification among those with the outcome. The first argument provides the probability distribution function (uniform, triangular, trapezoidal, logit-logistic, or logit-normal) and the second its parameters as a vector:
#' \enumerate{
#' \item Constant: constant value,
#' \item Uniform: min, max,
#' \item Triangular: lower limit, upper limit, mode,
#' \item Trapezoidal: min, lower mode, upper mode, max,
#' \item Logit-logistic: location, scale, lower bound shift, upper bound shift,
#' \item Logit-normal: location, scale, lower bound shift, upper bound shift.
#' }
#' @param seexp.parms List defining the sensitivity of exposure classification among those without the outcome.
#' @param spca.parms List defining the specificity of exposure classification among those with the outcome.
#' @param spexp.parms List defining the specifity of exposure classification among those without the outcome.
#' @param corr.se Correlation between case and non-case sensitivities.
#' @param corr.sp Correlation between case and non-case specificities.
#' @param discard A logical scalar. In case of negative adjusted count, should the draws be discarded? If set to FALSE, negative counts are set to zero.
#' @param alpha Significance level.
#' @param dec Number of decimals in the printout.
#' @param print A logical scalar. Should the results be printed?
#'
#' @return A list with elements:
#' \item{obs.data}{The analysed 2 x 2 table from the observed data.}
#' \item{obs.measures}{A table of observed incidence rate ratio with exact confidence interval.}
#' \item{adj.measures}{A table of corrected incidence rate ratios.}
#' \item{sim.df}{Data frame of random parameters and computed values.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#' @examples
#' set.seed(123)
#' # Exposure misclassification, non-differential
#' probsens.irr(matrix(c(2, 67232, 58, 10539000),
#' dimnames = list(c("GBS+", "Person-time"), c("HPV+", "HPV-")), ncol = 2),
#' reps = 20000,
#' seca.parms = list("trapezoidal", c(.4, .45, .55, .6)),
#' spca.parms = list("constant", 1))
#' @export
#' @importFrom stats binom.test median quantile qnorm runif
probsens.irr <- function(counts,
                         pt = NULL,
                         reps = 1000,
                         seca.parms = list(dist = c("constant", "uniform",
                                               "triangular", "trapezoidal",
                                               "logit-logistic", "logit-normal"),
                             parms = NULL),
                         seexp.parms = NULL,
                         spca.parms = list(dist = c("constant", "uniform",
                                               "triangular", "trapezoidal",
                                               "logit-logistic", "logit-normal"),
                             parms = NULL),
                         spexp.parms = NULL,
                         corr.se = NULL,
                         corr.sp = NULL,
                         discard = TRUE,
                         alpha = 0.05,
                         dec = 4,
                         print = TRUE){
    if(reps < 1)
        stop(paste("Invalid argument: reps = ", reps))
    
    if(is.null(seca.parms) | is.null(spca.parms))
        stop('At least one Se and one Sp should be provided through outcome parameters.')
    if(!is.list(seca.parms))
        stop('Sensitivity of exposure classification among those with the outcome should be a list.')
    else seca.parms <- seca.parms
    if(!is.null(corr.se) && (seca.parms[[1]] == "constant" | seexp.parms[[1]] == "constant"))
        stop('No correlated distributions with constant values.')
    if(!is.null(corr.sp) && (spca.parms[[1]] == "constant" | spexp.parms[[1]] == "constant"))
        stop('No correlated distributions with constant values.')
    if(seca.parms[[1]] == "constant" & length(seca.parms[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(seca.parms[[1]] == "uniform" & length(seca.parms[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(seca.parms[[1]] == "uniform" & seca.parms[[2]][1] >= seca.parms[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(seca.parms[[1]] == "triangular" & length(seca.parms[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(seca.parms[[1]] == "triangular" & ((seca.parms[[2]][1] > seca.parms[[2]][3]) |
                                        (seca.parms[[2]][2] < seca.parms[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(seca.parms[[1]] == "trapezoidal" & length(seca.parms[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(seca.parms[[1]] == "trapezoidal" & ((seca.parms[[2]][1] > seca.parms[[2]][2]) |
                                         (seca.parms[[2]][2] > seca.parms[[2]][3]) |
                                         (seca.parms[[2]][3] > seca.parms[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')
    if(seca.parms[[1]] == "logit-logistic" & (length(seca.parms[[2]]) < 2 | length(seca.parms[[2]]) == 3 | length(seca.parms[[2]]) > 4))
        stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(seca.parms[[1]] == "logit-logistic" & length(seca.parms[[2]]) == 4 &
       ((seca.parms[[2]][3] >= seca.parms[[2]][4]) | (!all(seca.parms[[2]][3:4] >= 0 & seca.parms[[2]][3:4] <= 1))))
        stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(seca.parms[[1]] == "logit-logistic" & length(seca.parms[[2]]) == 2)
        seca.parms <- list(seca.parms[[1]], c(seca.parms[[2]], c(0, 1)))
    if(seca.parms[[1]] == "logit-normal" & (length(seca.parms[[2]]) < 2 | length(seca.parms[[2]]) == 3 | length(seca.parms[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(seca.parms[[1]] == "logit-normal" & length(seca.parms[[2]]) == 4 &
       ((seca.parms[[2]][3] >= seca.parms[[2]][4]) | (!all(seca.parms[[2]][3:4] >= 0 & seca.parms[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(seca.parms[[1]] == "logit-normal" & length(seca.parms[[2]]) == 2)
        seca.parms <- list(seca.parms[[1]], c(seca.parms[[2]], c(0, 1)))
    if((seca.parms[[1]] == "constant" | seca.parms[[1]] == "uniform" | seca.parms[[1]] == "triangular" | seca.parms[[1]] == "trapezoidal") & !all(seca.parms[[2]] >= 0 & seca.parms[[2]] <= 1))
        stop('Sensitivity of exposure classification among those with the outcome should be between 0 and 1.')
    
    if(!is.null(seexp.parms) & !is.list(seexp.parms))
        stop('Sensitivity of exposure classification among those without the outcome should be a list.')
    else seexp.parms <- seexp.parms
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "constant" &
       length(seexp.parms[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "uniform" &
       length(seexp.parms[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "uniform" &&
       seexp.parms[[2]][1] >= seexp.parms[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "triangular" &
       length(seexp.parms[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "triangular" &&
       ((seexp.parms[[2]][1] > seexp.parms[[2]][3]) |
            (seexp.parms[[2]][2] < seexp.parms[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "trapezoidal" &
       length(seexp.parms[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "trapezoidal" &&
       ((seexp.parms[[2]][1] > seexp.parms[[2]][2]) |
            (seexp.parms[[2]][2] > seexp.parms[[2]][3]) |
                (seexp.parms[[2]][3] > seexp.parms[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "logit-logistic" & (length(seexp.parms[[2]]) < 2 | length(seexp.parms[[2]]) == 3 | length(seexp.parms[[2]]) > 4))
        stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "logit-logistic" & length(seexp.parms[[2]]) == 4 && ((seexp.parms[[2]][3] >= seexp.parms[[2]][4]) | (!all(seexp.parms[[2]][3:4] >= 0 & seexp.parms[[2]][3:4] <= 1))))
        stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "logit-logistic" & length(seexp.parms[[2]]) == 2)
        seexp.parms <- list(seexp.parms[[1]], c(seexp.parms[[2]], c(0, 1)))
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "logit-normal" & (length(seexp.parms[[2]]) < 2 | length(seexp.parms[[2]]) == 3 | length(seexp.parms[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "logit-normal" & length(seexp.parms[[2]]) == 4 && ((seexp.parms[[2]][3] >= seexp.parms[[2]][4]) | (!all(seexp.parms[[2]][3:4] >= 0 & seexp.parms[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(!is.null(seexp.parms) && seexp.parms[[1]] == "logit-normal" & length(seexp.parms[[2]]) == 2)
        seexp.parms <- list(seexp.parms[[1]], c(seexp.parms[[2]], c(0, 1)))
    if(!is.null(seexp.parms) && (seexp.parms[[1]] == "constant" |seexp.parms[[1]] == "uniform" | seexp.parms[[1]] == "triangular" | seexp.parms[[1]] == "trapezoidal") & !all(seexp.parms[[2]] >= 0 & seexp.parms[[2]] <= 1))
        stop('Sensitivity of exposure classification among those without the outcome should be between 0 and 1.')
    
    if(!is.list(spca.parms))
        stop('Specificity of exposure classification among those with the outcome should be a list.')
    else spca.parms <- spca.parms
    if(spca.parms[[1]] == "constant" & length(spca.parms[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(spca.parms[[1]] == "uniform" & length(spca.parms[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(spca.parms[[1]] == "uniform" & spca.parms[[2]][1] >= spca.parms[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(spca.parms[[1]] == "triangular" & length(spca.parms[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(spca.parms[[1]] == "triangular" & ((spca.parms[[2]][1] > spca.parms[[2]][3]) |
                                        (spca.parms[[2]][2] < spca.parms[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(spca.parms[[1]] == "trapezoidal" & length(spca.parms[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(spca.parms[[1]] == "trapezoidal" & ((spca.parms[[2]][1] > spca.parms[[2]][2]) |
                                         (spca.parms[[2]][2] > spca.parms[[2]][3]) |
                                         (spca.parms[[2]][3] > spca.parms[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')    
    if(spca.parms[[1]] == "logit-logistic" & (length(spca.parms[[2]]) < 2 | length(spca.parms[[2]]) == 3 | length(spca.parms[[2]]) > 4))
        stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(spca.parms[[1]] == "logit-logistic" & length(spca.parms[[2]]) == 4 &
       ((spca.parms[[2]][3] >= spca.parms[[2]][4]) | (!all(spca.parms[[2]][3:4] >= 0 & spca.parms[[2]][3:4] <= 1))))
        stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(spca.parms[[1]] == "logit-logistic" & length(spca.parms[[2]]) == 2)
        spca.parms <- list(spca.parms[[1]], c(spca.parms[[2]], c(0, 1)))
    if(spca.parms[[1]] == "logit-normal" & (length(spca.parms[[2]]) < 2 | length(spca.parms[[2]]) == 3 | length(spca.parms[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(spca.parms[[1]] == "logit-normal" & length(spca.parms[[2]]) == 4 &
       ((spca.parms[[2]][3] >= spca.parms[[2]][4]) | (!all(spca.parms[[2]][3:4] >= 0 & spca.parms[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(spca.parms[[1]] == "logit-normal" & length(spca.parms[[2]]) == 2)
        spca.parms <- list(spca.parms[[1]], c(spca.parms[[2]], c(0, 1)))
    if((spca.parms[[1]] == "constant" | spca.parms[[1]] == "uniform" | spca.parms[[1]] == "triangular" | spca.parms[[1]] == "trapezoidal") & !all(spca.parms[[2]] >= 0 & spca.parms[[2]] <= 1))
        stop('Specificity of exposure classification among those with the outcome should be between 0 and 1.')
    
    if(!is.null(spexp.parms) & !is.list(spexp.parms))
        stop('Specificity of exposure classification among those without the outcome should be a list.')
    else spexp.parms <- spexp.parms
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "constant" &
       length(spexp.parms[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "uniform" &
       length(spexp.parms[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "uniform" &&
       spexp.parms[[2]][1] >= spexp.parms[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "triangular" &
       length(spexp.parms[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "triangular" &&
       ((spexp.parms[[2]][1] > spexp.parms[[2]][3]) |
            (spexp.parms[[2]][2] < spexp.parms[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "trapezoidal" &
       length(spexp.parms[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "trapezoidal" &&
       ((spexp.parms[[2]][1] > spexp.parms[[2]][2]) |
            (spexp.parms[[2]][2] > spexp.parms[[2]][3]) |
                (spexp.parms[[2]][3] > spexp.parms[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')    
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "logit-logistic" & (length(spexp.parms[[2]]) < 2 | length(spexp.parms[[2]]) == 3 | length(spexp.parms[[2]]) > 4))
        stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "logit-logistic" & length(spexp.parms[[2]]) == 4 && ((spexp.parms[[2]][3] >= spexp.parms[[2]][4]) | (!all(spexp.parms[[2]][3:4] >= 0 & spexp.parms[[2]][3:4] <= 1))))
        stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(!is.null(seexp.parms) && spexp.parms[[1]] == "logit-logistic" & length(spexp.parms[[2]]) == 2)
        spexp.parms <- list(spexp.parms[[1]], c(spexp.parms[[2]], c(0, 1)))
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "logit-normal" & (length(spexp.parms[[2]]) < 2 | length(spexp.parms[[2]]) == 3 | length(spexp.parms[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "logit-normal" & length(spexp.parms[[2]]) == 4 && ((spexp.parms[[2]][3] >= spexp.parms[[2]][4]) | (!all(spexp.parms[[2]][3:4] >= 0 & spexp.parms[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(!is.null(spexp.parms) && spexp.parms[[1]] == "logit-normal" & length(spexp.parms[[2]]) == 2)
        spexp.parms <- list(spexp.parms[[1]], c(spexp.parms[[2]], c(0, 1)))
    if(!is.null(spexp.parms) && (spexp.parms[[1]] == "constant" | spexp.parms[[1]] == "uniform" | spexp.parms[[1]] == "triangular" | spexp.parms[[1]] == "trapezoidal") & !all(spexp.parms[[2]] >= 0 & spexp.parms[[2]] <= 1))
        stop('Specificity of exposure classification among those without the outcome should be between 0 and 1.')
    
    if(!is.null(seexp.parms) & (is.null(spca.parms) | is.null(spexp.parms) |
                                is.null(corr.se) | is.null(corr.sp)))
        stop('For non-differential misclassification type, have to provide Se and Sp for among those with and without the outcome as well as Se and Sp correlations.')

    if(!is.null(corr.se) && (corr.se == 0 | corr.se == 1))
        stop('Correlations should be > 0 and < 1.')
    if(!is.null(corr.sp) && (corr.sp == 0 | corr.sp == 1))
        stop('Correlations should be > 0 and < 1.')

    if(!is.null(pt) && inherits(counts, c("table", "matrix")))
        stop("pt argument should be NULL.")
    if(!inherits(counts, c("vector", "table", "matrix")))
        stop("counts argument should be a vector, a table, or a matrix.")
    if(is.null(pt) && inherits(counts, c("table", "matrix")))
        tab <- counts
    else tab <- rbind(counts, pt)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    draws <- matrix(NA, nrow = reps, ncol = 11)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A1", "B1", "C1", "D1",
                         "corr.IRR", "tot.IRR", "reps")
    corr.draws <- matrix(NA, nrow = reps, ncol = 10)

    seca <- c(reps, seca.parms[[2]])
    seexp <- c(reps, seexp.parms[[2]])
    spca <- c(reps, spca.parms[[2]])
    spexp <- c(reps, spexp.parms[[2]])

    obs.irr <- (a / c) / (b / d)
    lci.obs.irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1]) * c)
    uci.obs.irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2]) * c)

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
    
    if (is.null(seexp.parms) & !is.null(spca.parms) & is.null(spexp.parms) &
        is.null(corr.se) & is.null(corr.sp)) {
        if (seca.parms[[1]] == "constant") {
            draws[, 1] <- seca.parms[[2]]
        }
        if (seca.parms[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(seca))
            }
        if (seca.parms[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(seca))
            }
        if (seca.parms[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(seca))
            }
        if (seca.parms[[1]] == "logit-logistic") {
            draws[, 1] <- logitlog.dstr(seca)
            }
        if (seca.parms[[1]] == "logit-normal") {
            draws[, 1] <- logitnorm.dstr(seca)
            }
        draws[, 2] <- draws[, 1]
        if (spca.parms[[1]] == "constant") {
            draws[, 3] <- spca.parms[[2]]
        }
        if (spca.parms[[1]] == "uniform") {
            draws[, 3] <- do.call(runif, as.list(spca))
            }
        if (spca.parms[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::rtriangle, as.list(spca))
            }
        if (spca.parms[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(spca))
            }
        if (spca.parms[[1]] == "logit-logistic") {
            draws[, 3] <- logitlog.dstr(spca)
            }
        if (spca.parms[[1]] == "logit-normal") {
            draws[, 3] <- logitnorm.dstr(spca)
            }
        draws[, 4] <- draws[, 3]
    } else {
        corr.draws[, 1:6] <- apply(corr.draws[, 1:6],
                                   2,
                                   function(x) x = runif(reps))
        corr.draws[, 1:6] <- apply(corr.draws[, 1:6],
                                   2,
                                   function(x) log(x / (1 - x)))
        corr.draws[, 7] <- exp(sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 2]) /
            (1 + (exp(sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 2])))
        corr.draws[, 8] <- exp(sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 3]) /
            (1 + (exp(sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 3])))
        corr.draws[, 9] <- exp(sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 5]) /
            (1 + (exp(sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 5])))
        corr.draws[, 10] <- exp(sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 6]) /
            (1 + (exp(sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 6])))

    if (seca.parms[[1]] == "uniform") {
        draws[, 1] <- seca.parms[[2]][2] -
            (seca.parms[[2]][2] - seca.parms[[2]][1]) * corr.draws[, 7]
    }
    if (seca.parms[[1]] == "triangular") {
        draws[, 1] <- (corr.draws[, 7] *
            (seca.parms[[2]][2] - seca.parms[[2]][1]) + (seca.parms[[2]][1] + seca.parms[[2]][3])) / 2
        draws[, 1] <- ifelse(draws[, 1] < seca.parms[[2]][3],
                             seca.parms[[2]][1] + sqrt(abs((seca.parms[[2]][3] - seca.parms[[2]][1]) * (2 * draws[, 1] - seca.parms[[2]][1] - seca.parms[[2]][3]))),
                             draws[, 1])
        draws[, 1] <- ifelse(draws[, 1] > seca.parms[[2]][3],
                             seca.parms[[2]][2] - sqrt(abs(2 * (seca.parms[[2]][2] - seca.parms[[2]][3]) * (draws[, 1] - seca.parms[[2]][3]))),
                             draws[, 1])
    }
    if (seca.parms[[1]] == "trapezoidal") {
        draws[, 1] <- (corr.draws[, 7] *
            (seca.parms[[2]][4] + seca.parms[[2]][3] - seca.parms[[2]][1] - seca.parms[[2]][2]) + (seca.parms[[2]][1] + seca.parms[[2]][2])) / 2
        draws[, 1] <- ifelse(draws[, 1] < seca.parms[[2]][2],
                             seca.parms[[2]][1] + sqrt(abs((seca.parms[[2]][2] - seca.parms[[2]][1]) * (2 * draws[, 1] - seca.parms[[2]][1] - seca.parms[[2]][2]))),
                             draws[, 1])
        draws[, 1] <- ifelse(draws[, 1] > seca.parms[[2]][3],
                             seca.parms[[2]][4] - sqrt(abs(2 * (seca.parms[[2]][4] - seca.parms[[2]][3]) * (draws[, 1] - seca.parms[[2]][3]))),
                             draws[, 1])
    }
    if (seca.parms[[1]] == "logit-logistic") {
        seca.w <- seca.parms[[2]][1] + (seca.parms[[2]][2] * log(corr.draws[, 7] / (1 - corr.draws[, 7])))
        draws[, 1] <- seca.parms[[2]][3] + (seca.parms[[2]][4] - seca.parms[[2]][3]) * exp(seca.w) / (1 + exp(seca.w))
    }
    if (seca.parms[[1]] == "logit-normal") {
        seca.w <- seca.parms[[2]][1] + (seca.parms[[2]][2] * qnorm(corr.draws[, 7]))
        draws[, 1] <- seca.parms[[2]][3] + (seca.parms[[2]][4] - seca.parms[[2]][3]) * exp(seca.w) / (1 + exp(seca.w))
    }
    if (seexp.parms[[1]] == "uniform") {
        draws[, 2] <- seexp.parms[[2]][2] -
            (seexp.parms[[2]][2] - seexp.parms[[2]][1]) * corr.draws[, 8]
    }
    if (seexp.parms[[1]] == "triangular") {
        draws[, 2] <- (corr.draws[, 8] *
                           (seexp.parms[[2]][2] - seexp.parms[[2]][1]) + (seexp.parms[[2]][1] + seexp.parms[[2]][3])) / 2
        draws[, 2] <- ifelse(draws[, 2] < seexp.parms[[2]][3],
                             seexp.parms[[2]][1] + sqrt(abs((seexp.parms[[2]][3] - seexp.parms[[2]][1]) * (2 * draws[, 2] - seexp.parms[[2]][1] - seexp.parms[[2]][3]))),
                             draws[, 2])
        draws[, 2] <- ifelse(draws[, 2] > seexp.parms[[2]][3],
                             seexp.parms[[2]][2] - sqrt(abs(2 * (seexp.parms[[2]][2] - seexp.parms[[2]][3]) * (draws[, 2] - seexp.parms[[2]][3]))),
                             draws[, 2])
    }
    if (seexp.parms[[1]] == "trapezoidal") {
        draws[, 2] <- (corr.draws[, 8] *
                           (seexp.parms[[2]][4] + seexp.parms[[2]][3] - seexp.parms[[2]][1] - seexp.parms[[2]][2]) + (seexp.parms[[2]][1] + seexp.parms[[2]][2])) / 2
        draws[, 2] <- ifelse(draws[, 2] < seexp.parms[[2]][2],
                             seexp.parms[[2]][1] + sqrt(abs((seexp.parms[[2]][2] - seexp.parms[[2]][1]) * (2 * draws[, 2] - seexp.parms[[2]][1] - seexp.parms[[2]][2]))),
                             draws[, 2])
        draws[, 2] <- ifelse(draws[, 2] > seexp.parms[[2]][3],
                             seexp.parms[[2]][4] - sqrt(abs(2 * (seexp.parms[[2]][4] - seexp.parms[[2]][3]) * (draws[, 2] - seexp.parms[[2]][3]))),
                             draws[, 2])
    }
    if (seexp.parms[[1]] == "logit-logistic") {
        seexp.w <- seexp.parms[[2]][1] + (seexp.parms[[2]][2] * log(corr.draws[, 8] / (1 - corr.draws[, 8])))
        draws[, 2] <- seexp.parms[[2]][3] + (seexp.parms[[2]][4] - seexp.parms[[2]][3]) * exp(seexp.w) / (1 + exp(seexp.w))
    }
    if (seexp.parms[[1]] == "logit-normal") {
        seexp.w <- seexp.parms[[2]][1] + (seexp.parms[[2]][2] * qnorm(corr.draws[, 8]))
        draws[, 2] <- seexp.parms[[2]][3] + (seexp.parms[[2]][4] - seexp.parms[[2]][3]) * exp(seexp.w) / (1 + exp(seexp.w))
    }
    if (spca.parms[[1]] == "uniform") {
        draws[, 3] <- spca.parms[[2]][2] -
            (spca.parms[[2]][2] - spca.parms[[2]][1]) * corr.draws[, 9]
    }
    if (spca.parms[[1]] == "triangular") {
        draws[, 3] <- (corr.draws[, 9] *
                           (spca.parms[[2]][2] - spca.parms[[2]][1]) + (spca.parms[[2]][1] + spca.parms[[2]][3])) / 2
        draws[, 3] <- ifelse(draws[, 3] < spca.parms[[2]][3],
                             spca.parms[[2]][1] + sqrt(abs((spca.parms[[2]][3] - spca.parms[[2]][1]) * (2 * draws[, 3] - spca.parms[[2]][1] - spca.parms[[2]][3]))),
                             draws[, 3])
        draws[, 3] <- ifelse(draws[, 3] > spca.parms[[2]][3],
                             spca.parms[[2]][2] - sqrt(abs(2 * (spca.parms[[2]][2] - spca.parms[[2]][3]) * (draws[, 3] - spca.parms[[2]][3]))),
                             draws[, 3])
    }
    if (spca.parms[[1]] == "trapezoidal") {
        draws[, 3] <- (corr.draws[, 9] *
                           (spca.parms[[2]][4] + spca.parms[[2]][3] - spca.parms[[2]][1] - spca.parms[[2]][2]) + (spca.parms[[2]][1] + spca.parms[[2]][2])) / 2
        draws[, 3] <- ifelse(draws[, 3] < spca.parms[[2]][2],
                             spca.parms[[2]][1] + sqrt(abs((spca.parms[[2]][2] - spca.parms[[2]][1]) * (2 * draws[, 3] - spca.parms[[2]][1] - spca.parms[[2]][2]))),
                             draws[, 3])
        draws[, 3] <- ifelse(draws[, 3] > spca.parms[[2]][3],
                             spca.parms[[2]][4] - sqrt(abs(2 * (spca.parms[[2]][4] - spca.parms[[2]][3]) * (draws[, 3] - spca.parms[[2]][3]))),
                             draws[, 3])
    }
    if (spca.parms[[1]] == "logit-logistic") {
        spca.w <- spca.parms[[2]][1] + (spca.parms[[2]][2] * log(corr.draws[, 9] / (1 - corr.draws[, 9])))
        draws[, 3] <- spca.parms[[2]][3] + (spca.parms[[2]][4] - spca.parms[[2]][3]) * exp(spca.w) / (1 + exp(spca.w))
    }
    if (spca.parms[[1]] == "logit-normal") {
        spca.w <- spca.parms[[2]][1] + (spca.parms[[2]][2] * qnorm(corr.draws[, 9]))
        draws[, 3] <- spca.parms[[2]][3] + (spca.parms[[2]][4] - spca.parms[[2]][3]) * exp(spca.w) / (1 + exp(spca.w))
    }
    if (spexp.parms[[1]] == "uniform") {
        draws[, 4] <- spexp.parms[[2]][2] -
            (spexp.parms[[2]][2] - spexp.parms[[2]][1]) * corr.draws[, 10]
    }
    if (spexp.parms[[1]] == "triangular") {
        draws[, 4] <- (corr.draws[, 10] *
                           (spexp.parms[[2]][2] - spexp.parms[[2]][1]) + (spexp.parms[[2]][1] + spexp.parms[[2]][3])) / 2
        draws[, 4] <- ifelse(draws[, 4] < spexp.parms[[2]][3],
                             spexp.parms[[2]][1] + sqrt(abs((spexp.parms[[2]][3] - spexp.parms[[2]][1]) * (2 * draws[, 4] - spexp.parms[[2]][1] - spexp.parms[[2]][3]))),
                             draws[, 4])
        draws[, 4] <- ifelse(draws[, 4] > spexp.parms[[2]][3],
                             spexp.parms[[2]][2] - sqrt(abs(2 * (spexp.parms[[2]][2] - spexp.parms[[2]][3]) * (draws[, 4] - spexp.parms[[2]][3]))),
                             draws[, 4])
    }
    if (spexp.parms[[1]] == "trapezoidal") {
        draws[, 4] <- (corr.draws[, 10] *
                           (spexp.parms[[2]][4] + spexp.parms[[2]][3] - spexp.parms[[2]][1] - spexp.parms[[2]][2]) + (spexp.parms[[2]][1] + spexp.parms[[2]][2])) / 2
        draws[, 4] <- ifelse(draws[, 4] < spexp.parms[[2]][2],
                             spexp.parms[[2]][1] + sqrt(abs((spexp.parms[[2]][2] - spexp.parms[[2]][1]) * (2 * draws[, 4] - spexp.parms[[2]][1] - spexp.parms[[2]][2]))),
                             draws[, 4])
        draws[, 4] <- ifelse(draws[, 4] > spexp.parms[[2]][3],
                             spexp.parms[[2]][4] - sqrt(abs(2 * (spexp.parms[[2]][4] - spexp.parms[[2]][3]) * (draws[, 4] - spexp.parms[[2]][3]))),
                             draws[, 4])
    }
    if (spexp.parms[[1]] == "logit-logistic") {
        spexp.w <- spexp.parms[[2]][1] + (spexp.parms[[2]][2] * log(corr.draws[, 10] / (1 - corr.draws[, 10])))
        draws[, 4] <- spexp.parms[[2]][3] + (spexp.parms[[2]][4] - spexp.parms[[2]][3]) * exp(spexp.w) / (1 + exp(spexp.w))
    }
    if (spexp.parms[[1]] == "logit-normal") {
        spexp.w <- spexp.parms[[2]][1] + (spexp.parms[[2]][2] * qnorm(corr.draws[, 10]))
        draws[, 4] <- spexp.parms[[2]][3] + (spexp.parms[[2]][4] - spexp.parms[[2]][3]) * exp(spexp.w) / (1 + exp(spexp.w))
    }
    }
    
    draws[, 11] <- runif(reps)

    draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) /
        (draws[, 1] - (1 - draws[, 3]))
    draws[, 6] <- (a + b) - draws[, 5]
    draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) /
        (draws[, 2] - (1 - draws[, 4]))
    draws[, 8] <- (c + d) - draws[, 7]
    
    draws[, 9] <- (draws[, 5]/(draws[, 5] + draws[, 7])) /
        (draws[, 6]/(draws[, 6] + draws[, 8]))
    
    draws[, 9] <- ifelse(draws[, 5] < 1 |
                             draws[, 6] < 1 |
                                 draws[, 7] < 1 |
                                     draws[, 8] < 1, NA, draws[, 9])
    draws[, 9] <- ifelse(draws[, 2] < (c / (c + d)) |
                             draws[, 4] < (d / (c + d)) |
                                 draws[, 1] < (a / (a + b)) |
                                     draws[, 3] < (b / (a + b)), NA, draws[, 9])
    if(all(is.na(draws[, 9])))
        warning('Prior Se/Sp distributions lead to all negative adjusted counts.')
    if (discard) {
        if(sum(is.na(draws[, 9])) > 0)
            message('Chosen prior Se/Sp distributions lead to ',
                    sum(is.na(draws[, 9])),
                    ' negative adjusted counts which were discarded.')
    }
    else {
        if(sum(is.na(draws[, 9])) > 0) {
            message('Chosen prior Se/Sp distributions lead to ',
                    sum(is.na(draws[, 9])),
                    ' negative adjusted counts which were set to zero.')
        }
        draws[, 9] <- ifelse(is.na(draws[, 9]), 0, draws[, 9])
    }

    draws[, 10] <- exp(log(draws[, 9]) -
                           qnorm(draws[, 11]) *
                               ((log(uci.obs.irr) - log(lci.obs.irr)) /
                                    (qnorm(.975) * 2)))
    
    corr.irr <- c(median(draws[, 9], na.rm = TRUE),
                 quantile(draws[, 9], probs = .025, na.rm = TRUE),
                 quantile(draws[, 9], probs = .975, na.rm = TRUE))
    tot.irr <- c(median(draws[, 10], na.rm = TRUE),
                quantile(draws[, 10], probs = .025, na.rm = TRUE),
                quantile(draws[, 10], probs = .975, na.rm = TRUE))

    if (is.null(rownames(tab)))
        rownames(tab) <- c("Cases", "Person-time")
    if (is.null(colnames(tab)))
        colnames(tab) <- c("Exposed", "Unexposed")
    if (print) {
        cat("Observed Data:",
            "\n-------------------------",
            "\n\n")
        print(round(tab, dec))
        cat("\n")
    }
    rmat <- data.frame(obs.irr, lci.obs.irr, uci.obs.irr)
    rownames(rmat) <- " Observed Incidence Rate ratio:"
    colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                       sep = ""), "interval")
    if (print) {
        cat("Observed Measures:",
            "\n-----------------------------------------------------\n\n")
        print(round(rmat, dec))
        cat("\n")
    }
    rmatc <- rbind(corr.irr, tot.irr)
    rownames(rmatc) <- c("           Incidence Rate Ratio -- systematic error:",
                         "Incidence Rate Ratio -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    if (print) {
        print(round(rmatc, dec))
        cat("\nBias Parameters:",
            "\n----------------\n\n")
        cat("   Se|Cases:", seca.parms[[1]], "(", seca.parms[[2]], ")",
            "\n   Sp|Cases:", spca.parms[[1]], "(", spca.parms[[2]], ")",
            "\nSe|No-cases:", seexp.parms[[1]], "(", seexp.parms[[2]], ")",
            "\nSp|No-cases:", spexp.parms[[1]], "(", spexp.parms[[2]], ")",
            "\nDiscard negative adjusted counts:", discard,
            "\n")
    }
    invisible(list(obs.data = tab,
                   obs.measures = rmat, 
                   adj.measures = rmatc,
                   sim.df = as.data.frame(draws[, -11])))
}
