#' Probabilistic sensitivity analysis.
#'
#' Probabilistic sensitivity analysis to correct for exposure misclassification or
#' outcome misclassification and random error.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param type Choice of correction for exposure or outcome misclassification.
#' @param reps Number of replications to run.
#' @param seca.parms List defining the sensitivity of exposure classification among those with the outcome. The first argument provides the probability distribution function (constant, uniform, triangular, trapezoidal, logit-logistic, or logit-normal) and the second its parameters as a vector:
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
#' \item{obs.measures}{A table of observed relative risk and odds ratio with confidence intervals.}
#' \item{adj.measures}{A table of corrected relative risks and odds ratios.}
#' \item{sim.df}{Data frame of random parameters and computed values.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#' @examples
#' # The data for this example come from:
#' # Greenland S., Salvan A., Wegman D.H., Hallock M.F., Smith T.J.
#' # A case-control study of cancer mortality at a transformer-assembly facility.
#' # Int Arch Occup Environ Health 1994; 66(1):49-54.
#' set.seed(123)
#' # Exposure misclassification, non-differential
#' probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 20000,
#' seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
#' spca.parms = list("trapezoidal", c(.75, .85, .95, 1)))
#' # Exposure misclassification, differential
#' probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 20000,
#' seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
#' seexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
#' spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
#' spexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
#' corr.se = .8,
#' corr.sp = .8)
#' # Disease misclassification
#' probsens(matrix(c(173, 602, 134, 663),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "outcome",
#' reps = 20000,
#' seca.parms = list("uniform", c(.8, 1)),
#' spca.parms = list("uniform", c(.8, 1)))
#' @export
#' @importFrom stats median qnorm quantile runif
probsens <- function(case,
                     exposed,
                     type = c("exposure", "outcome"),
                     reps = 1000,
                     seca.parms = list(dist = c("constant", "uniform", "triangular",
                                           "trapezoidal", "logit-logistic",
                                                "logit-normal"),
                                       parms = NULL),
                     seexp.parms = NULL,
                     spca.parms = list(dist = c("constant", "uniform", "triangular",
                                           "trapezoidal", "logit-logistic",
                                                "logit-normal"),
                                       parms = NULL),
                     spexp.parms = NULL,
                     corr.se = NULL,
                     corr.sp = NULL,
                     discard = TRUE,
                     alpha = 0.05,
                     dec = 4,
                     print = TRUE){
    if(reps < 1)
        stop(paste("Invalid argument: reps =", reps))
    
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
    if(!is.null(seexp.parms) && (seexp.parms[[1]] == "constant" | seexp.parms[[1]] == "uniform" | seexp.parms[[1]] == "triangular" | seexp.parms[[1]] == "trapezoidal") & !all(seexp.parms[[2]] >= 0 & seexp.parms[[2]] <= 1))
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
    
    if(inherits(case, c("table", "matrix")))
        tab <- case
    else tab <- table(case, exposed)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    draws <- matrix(NA, nrow = reps, ncol = 13)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A1", "B1", "C1", "D1",
                         "corr.RR", "corr.OR",
                         "reps",
                         "tot.RR", "tot.OR")
    corr.draws <- matrix(NA, nrow = reps, ncol = 10)

    seca <- c(reps, seca.parms[[2]])
    seexp <- c(reps, seexp.parms[[2]])
    spca <- c(reps, spca.parms[[2]])
    spexp <- c(reps, spexp.parms[[2]])

    obs.rr <- (a/(a + c)) / (b/(b + d))
    se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
    lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
    uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)

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
            draws[, 2] <- spca.parms[[2]]
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

    type <- match.arg(type)
    if (type == "exposure") {
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) /
            (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (a + b) - draws[, 5]
        draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) /
            (draws[, 2] - (1 - draws[, 4]))
        draws[, 8] <- (c + d) - draws[, 7]

        draws[, 9] <- (draws[, 5]/(draws[, 5] + draws[, 7])) /
            (draws[, 6]/(draws[, 6] + draws[, 8]))
        draws[, 10] <- (draws[, 5]/draws[, 7]) / (draws[, 6]/draws[, 8])

        draws[, 9] <- ifelse(draws[, 5] < 1 |
                                 draws[, 6] < 1 |
                                     draws[, 7] < 1 |
                                         draws[, 8] < 1, NA, draws[, 9])
        draws[, 10] <- ifelse(draws[, 5] < 1 |
                                 draws[, 6] < 1 |
                                     draws[, 7] < 1 |
                                         draws[, 8] < 1, NA, draws[, 10])
        draws[, 9] <- ifelse(draws[, 2] < (c / (c + d)) |
                                 draws[, 4] < (d / (c + d)) |
                                     draws[, 1] < (a / (a + b)) |
                                         draws[, 3] < (b / (a + b)), NA, draws[, 9])
        draws[, 10] <- ifelse(draws[, 2] < (c / (c + d)) |
                                  draws[, 4] < (d / (c + d)) |
                                      draws[, 1] < (a / (a + b)) |
                                          draws[, 3] < (b / (a + b)), NA, draws[, 10])
        if(all(is.na(draws[, 9])) | all(is.na(draws[, 10])))
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
            draws[, 10] <- ifelse(is.na(draws[, 10]), 0, draws[, 10])
        }

        draws[, 12] <- exp(log(draws[, 9]) -
                               qnorm(draws[, 11]) *
                                         ((log(uci.obs.rr) - log(lci.obs.rr)) /
                                              (qnorm(.975) * 2)))
        draws[, 13] <- exp(log(draws[, 10]) -
                               qnorm(draws[, 11]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))

        corr.rr <- c(median(draws[, 9], na.rm = TRUE),
                     quantile(draws[, 9], probs = .025, na.rm = TRUE),
                     quantile(draws[, 9], probs = .975, na.rm = TRUE))
        corr.or <- c(median(draws[, 10], na.rm = TRUE),
                     quantile(draws[, 10], probs = .025, na.rm = TRUE),
                     quantile(draws[, 10], probs = .975, na.rm = TRUE))
        tot.rr <- c(median(draws[, 12], na.rm = TRUE),
                     quantile(draws[, 12], probs = .025, na.rm = TRUE),
                     quantile(draws[, 12], probs = .975, na.rm = TRUE))
        tot.or <- c(median(draws[, 13], na.rm = TRUE),
                     quantile(draws[, 13], probs = .025, na.rm = TRUE),
                     quantile(draws[, 13], probs = .975, na.rm = TRUE))        

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (print) {
            cat("Observed Data:",
                "\n--------------", 
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
            print(round(tab, dec))
            cat("\n")
            }
        rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr),
                      c(obs.or, lci.obs.or, uci.obs.or))
        rownames(rmat) <- c(" Observed Relative Risk:", "    Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                               sep = ""), "interval")
        if (print) {
            cat("Observed Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
            print(round(rmat, dec))
            cat("\n")
            }
        rmatc <- rbind(corr.rr, corr.or, tot.rr, tot.or)
        rownames(rmatc) <- c("           Relative Risk -- systematic error:",
                             "              Odds Ratio -- systematic error:",
                             "Relative Risk -- systematic and random error:",
                             "   Odds Ratio -- systematic and random error:")
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
        }

    if (type == "outcome") {
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + c)) /
            (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (b - (1 - draws[, 4]) * (b + d)) /
            (draws[, 2] - (1 - draws[, 4]))
        draws[, 7] <- (a + c) - draws[, 5]
        draws[, 8] <- (b + d) - draws[, 6]

        draws[, 9] <- (draws[, 5]/(draws[, 5] + draws[, 7])) /
            (draws[, 6]/(draws[, 6] + draws[, 8]))
        draws[, 10] <- (draws[, 5]/draws[, 7]) / (draws[, 6]/draws[, 8])

        draws[, 9] <- ifelse(draws[, 5] < 1 |
                               draws[, 6] < 1 |
                                 draws[, 7] < 1 |
                                   draws[, 8] < 1, NA, draws[, 9])
        draws[, 10] <- ifelse(draws[, 5] < 1 |
                               draws[, 6] < 1 |
                                 draws[, 7] < 1 |
                                   draws[, 8] < 1, NA, draws[, 10])
        draws[, 9] <- ifelse(draws[, 2] < (c / (c + d)) |
                             draws[, 4] < (d / (c + d)) |
                             draws[, 1] < (a / (a + b)) |
                             draws[, 3] < (b / (a + b)), NA, draws[, 9])
        draws[, 10] <- ifelse(draws[, 2] < (c / (c + d)) |
                              draws[, 4] < (d / (c + d)) |
                              draws[, 1] < (a / (a + b)) |
                              draws[, 3] < (b / (a + b)), NA, draws[, 10])
        if(all(is.na(draws[, 9])) | all(is.na(draws[, 10])))
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
            draws[, 10] <- ifelse(is.na(draws[, 10]), 0, draws[, 10])
        }

        draws[, 12] <- exp(log(draws[, 9]) -
                               qnorm(draws[, 11]) *
                                         ((log(uci.obs.rr) - log(lci.obs.rr)) /
                                              (qnorm(.975) * 2)))
        draws[, 13] <- exp(log(draws[, 10]) -
                               qnorm(draws[, 11]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))

        corr.rr <- c(median(draws[, 9], na.rm = TRUE),
                     quantile(draws[, 9], probs = .025, na.rm = TRUE),
                     quantile(draws[, 9], probs = .975, na.rm = TRUE))
        corr.or <- c(median(draws[, 10], na.rm = TRUE),
                     quantile(draws[, 10], probs = .025, na.rm = TRUE),
                     quantile(draws[, 10], probs = .975, na.rm = TRUE))
        tot.rr <- c(median(draws[, 12], na.rm = TRUE),
                     quantile(draws[, 12], probs = .025, na.rm = TRUE),
                     quantile(draws[, 12], probs = .975, na.rm = TRUE))
        tot.or <- c(median(draws[, 13], na.rm = TRUE),
                     quantile(draws[, 13], probs = .025, na.rm = TRUE),
                     quantile(draws[, 13], probs = .975, na.rm = TRUE))        

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (print) {
            cat("Observed Data:",
                "\n--------------", 
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
            print(round(tab, dec))
            cat("\n")
            }
        rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr),
                      c(obs.or, lci.obs.or, uci.obs.or))
        rownames(rmat) <- c(" Observed Relative Risk:", "    Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                               sep = ""), "interval")
        if (print) {
            cat("Observed Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
            print(round(rmat, dec))
            cat("\n")
            }
        rmatc <- rbind(corr.rr, corr.or, tot.rr, tot.or)
        rownames(rmatc) <- c("           Relative Risk -- systematic error:",
                             "              Odds Ratio -- systematic error:",
                             "Relative Risk -- systematic and random error:",
                             "   Odds Ratio -- systematic and random error:")
        colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
        if (print) {
            print(round(rmatc, dec))
            cat("\nBias Parameters:",
                "\n----------------\n\n")
            cat("    Se|Exposed:", seca.parms[[1]], "(", seca.parms[[2]], ")",
                "\n    Sp|Exposed:", spca.parms[[1]], "(", spca.parms[[2]], ")",
                "\nSe|Non-exposed:", seexp.parms[[1]], "(", seexp.parms[[2]], ")",
                "\nSp|Non-exposed:", spexp.parms[[1]], "(", spexp.parms[[2]], ")",
                "\nDiscard negative adjusted counts:", discard,
                "\n")
            }
        }
    invisible(list(obs.data = tab,
                   obs.measures = rmat, 
                   adj.measures = rmatc,
                   sim.df = as.data.frame(draws[, -11])))
}
