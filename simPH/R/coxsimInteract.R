#' Simulate quantities of interest for linear multiplicative interactions */
#' from Cox Proportional Hazards models
#'
#' \code{coxsimInteract} simulates quantities of interest for linear
#' multiplicative interactions using multivariate normal distributions.
#' These can be plotted with \code{\link{simGG}}.
#' @param obj a \code{\link{coxph}} class fitted model object with a linear
#' multiplicative interaction.
#' @param b1 character string of the first constitutive variable's name.
#' Note \code{b1} and \code{b2} must be entered in the order in which they are
#' entered into the \code{coxph} model.
#' @param b2 character string of the second constitutive variable's name.
#' @param qi quantity of interest to simulate. Values can be
#' \code{"Marginal Effect"}, \code{"First Difference"}, \code{"Hazard Ratio"},
#' and \code{"Hazard Rate"}. The default is \code{qi = "Hazard Ratio"}.
#' If \code{qi = "Hazard Rate"} and the \code{coxph} model has strata, then
#' hazard rates for each strata will also be calculated.
#' @param X1 numeric vector of fitted values of \code{b1} to simulate for.
#' If \code{qi = "Marginal Effect"} then only \code{X2} can be set. If you want
#' to plot the results, \code{X1} should have more than one value.
#' @param X2 numeric vector of fitted values of \code{b2} to simulate for.
#' @param means logical, whether or not to use the mean values to fit the
#' hazard rate for covaraiates other than \code{b1} \code{b2} and \code{b1*b2}.
#' Note: it does not currently support models that include polynomials created
#' by \code{\link{I}}. Note: EXPERIMENTAL. \code{lines} are not currently
#' supported in \code{\link{simGG}} if \code{means = TRUE}.
#' @param expMarg logical. Whether or not to exponentiate the marginal effect.
#' @param nsim the number of simulations to run per value of X. Default is
#' \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is
#' \code{ci = 0.95}, i.e. keep the middle 95 percent. If \code{spin = TRUE}
#' then \code{ci} is the confidence level of the shortest probability interval.
#' Any value from 0 through 1 may be used.
#' @param spin logical, whether or not to keep only the shortest probability
#' interval rather than the middle simulations. Currently not supported for
#' hazard rates.
#' @param extremesDrop logical whether or not to drop simulated quantity of
#' interest values that are \code{Inf}, \code{NA}, \code{NaN} and
#' \eqn{> 1000000} for \code{spin = FALSE} or \eqn{> 800} for \code{spin = TRUE}.
#' These values are difficult to plot \code{\link{simGG}} and may prevent
#' \code{spin} from finding the central interval.
#'
#' @details Simulates marginal effects, first differences, hazard ratios, and
#' hazard rates for linear multiplicative interactions.
#' Marginal effects are calculated as in Brambor et al. (2006) with the
#' addition that we take the exponent, so that it resembles a hazard ratio.
#' You can choose not to take the exponent by setting the argument
#' \code{expMarg = FALSE}. For an interaction between variables \eqn{X} and
#' \eqn{Z} the marginal effect for \eqn{X} is:
#' \deqn{ME_{X} = e^(\beta_{X} + \beta_{XZ}Z)}{ME[X] = exp(\beta[X] +
#' \beta[XZ]Z)}
#'
#' Note that for First Differences the comparison is not between two values of
#' the same variable but two values of the constitute variable and 0 for the
#' two variables.
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ lethal*prevgenx,
#'             data = CarpenterFdaData)
#'
#' # Simulate Marginal Effect of lethal for multiple
#' # values of prevgenx
#' Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx",
#'                        X2 = seq(2, 115, by = 5), spin = TRUE)
#'
#' \dontrun{
#' # Change the order of the covariates to make a more easily
#' # interpretable relative hazard graph.
#' M2 <- coxph(Surv(acttime, censor) ~ prevgenx*lethal +
#'              orphdum, data = CarpenterFdaData)
#'
#' # Simulate Hazard Ratio of lethal for multiple
#' # values of prevgenx
#' Sim2 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
#'                     X1 = seq(2, 115, by = 2),
#'                     X2 = c(0, 1),
#'                     qi = "Hazard Ratio", ci = 0.9)
#'
#' # Simulate First Difference
#' Sim3 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
#'                        X1 = seq(2, 115, by = 2),
#'                        X2 = c(0, 1),
#'                        qi = "First Difference", spin = TRUE)
#'
#' # Simulate Hazard Rate
#' Sim4 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
#'                        X1 = 90, X2 = 1, qi = "Hazard Rate",
#'                        means = TRUE)
#' }
#'
#' # Example with a categorical variable
#' # Download data
#' library(rio)
#' hmohiv <- import("http://www.ats.ucla.edu/stat/r/examples/asa/hmohiv.csv")
#'
#' # Create category lables
#' hmohiv$drug <- factor(hmohiv$drug, labels = c('not treated', 'treated'))
#'
#' M3 <- coxph(Surv(time,censor) ~ drug*age, data = hmohiv)
#'
#' # Note: Use relevant coefficient name as shown in model summary, e.g.
#' # 'drugtreated'.
#' Sim5 <- coxsimInteract(M3, b1 = "drugtreated", b2 = 'age', X2 = 20:54)
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006.
#' ''Understanding Interaction Models: Improving Empirical Analyses.''
#' Political Analysis 14(1): 63-82.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. ''Making the Most of
#' Statistical Analyses: Improving Interpretation and Presentation.'' American
#' Journal of Political Science 44(2): 347-61.
#'
#' Liu, Ying, Andrew Gelman, and Tian Zheng. 2013. ''Simulation-Efficient
#' Shortest Probability Intervals.'' Arvix.
#' \url{http://arxiv.org/pdf/1302.2142v1.pdf}.
#'
#' @seealso \code{\link{simGG}}, \code{\link{survival}}, \code{\link{strata}},
#' and \code{\link{coxph}},
#' @return a \code{siminteract} class object
#' @import data.table
#' @importFrom stats vcov model.frame
#' @importFrom survival basehaz
#' @importFrom MASS mvrnorm
#' @export

coxsimInteract <- function(obj, b1, b2, qi = "Marginal Effect", X1 = NULL,
                           X2 = NULL, means = FALSE, expMarg = TRUE,
                           nsim = 1000, ci = 0.95, spin = FALSE,
                           extremesDrop = TRUE)
{
    HRValue <- strata <- QI <- SimID <- time <- NULL
    if (qi != "Hazard Rate" & isTRUE(means)) {
        stop("means can only be TRUE when qi = 'Hazard Rate'.", call. = FALSE)
    }
    # Ensure that qi is valid
    qiOpts <- c("Marginal Effect", "First Difference", "Hazard Ratio",
                "Hazard Rate")
    TestqiOpts <- qi %in% qiOpts
    if (!isTRUE(TestqiOpts)) {
        stop("Invalid qi type. qi must be 'Marginal Effect', 'First Difference', 'Hazard Ratio', or 'Hazard Rate'.",
            call. = FALSE)
    }
    MeansMessage <- NULL
    if (isTRUE(means) & length(obj$coefficients) == 3) {
        means <- FALSE
        MeansMessage <- FALSE
        message("Note: means reset to FALSE. The model only includes the interaction variables.")
    } else if (isTRUE(means) & length(obj$coefficients) > 3) {
        MeansMessage <- TRUE
    }

    # Create simulation ID variable
    SimID <- 1:nsim

    # Parameter estimates & Variance/Covariance matrix
    Coef <- matrix(obj$coefficients)
    VC <- vcov(obj)

    # Draw covariate estimates from the multivariate normal distribution
    Drawn <- mvrnorm(n = nsim, mu = Coef, Sigma = VC)
    DrawnDF <- data.frame(Drawn)
    dfn <- names(DrawnDF)

    bs <- c(b1, b2)
    bpos <- match(bs, dfn)
    binter <- paste0(bs[[1]], ".", bs[[2]])
    binter <- match(binter, dfn)
    NamesInt <- c(bpos, binter)

    # If all values aren't set for calculating the hazard rate
    if (!isTRUE(means)) {

        # Subset data frame to only include interaction constitutive terms and
        Simb <- data.frame(SimID, Drawn[, NamesInt])

        # Find quantity of interest
        if (qi == "Marginal Effect") {
            if (!is.null(X1)) {
                stop("For Marginal Effects only X2 should be specified.",
                    call. = FALSE)
            } else{
                X2df <- data.frame(X2)
                names(X2df) <- c("X2")
                Simb <- merge(Simb, X2df)
        if (isTRUE(expMarg)) {
                  Simb$QI <- exp(Simb[, 2] + (Simb[, 4] * Simb[, 5]))
        }
        else if (!isTRUE(expMarg)) {
          Simb$QI <- Simb[, 2] + (Simb[, 4] * Simb[, 5])
        }
            }
        }
        else if (qi == "First Difference") {
          if (is.null(X1) | is.null(X2)) {
            stop("For First Differences both X1 and X2 should be specified.",
                call. = FALSE)
          } else{
            Xs <- merge(X1, X2)
            names(Xs) <- c("X1", "X2")
            Xs$Comparison <- paste0(Xs[, 1], ", ", Xs[, 2])
            Simb <- merge(Simb, Xs)
            Simb$QI <- (exp((Simb$X1 * Simb[, 2]) + (Simb$X2 * Simb[, 3]) +
                       (Simb$X1 * Simb$X2 * Simb[, 4]) - 1) * 100)
          }
        }
        else if (qi == "Hazard Ratio") {
          if (is.null(X1) | is.null(X2)) {
            stop("For Hazard Ratios both X1 and X2 should be specified.",
                call. = FALSE)
          } else {
            Xs <- merge(X1, X2)
            names(Xs) <- c("X1", "X2")
            Xs$Comparison <- paste0(Xs[, 1], ", ", Xs[, 2])
            Simb <- merge(Simb, Xs)
            Simb$QI <- (exp((Simb$X1 * Simb[, 2]) + (Simb$X2 * Simb[, 3]) +
                       (Simb$X1 * Simb$X2 * Simb[, 4])))
          }
        }
        else if (qi == "Hazard Rate") {
            if (is.null(X1) | is.null(X2)) {
                stop("For Hazard Rates, both X1 and X2 should be specified.",
                     call. = FALSE)
            }
            if (isTRUE(MeansMessage)) {
                message("All variables' values other than b1, b2, and b1*b2 are fitted at 0.")
            }
            Xs <- data.frame(X1, X2)
            Xs$HRValue <- paste0(Xs$X1, ", ", Xs$X2)

            Simb <- merge(Simb, Xs)
            Simb$HR <- exp((Simb$X1 * Simb[, 2]) + (Simb$X2 * Simb[, 3]) +
                          (Simb$X1 * Simb$X2 * Simb[, 4]))

            bfit <- basehaz(obj)
            bfit$FakeID <- 1
            Simb$FakeID <- 1
            bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
            SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
            Simb <- SimbDT[bfitDT, allow.cartesian = TRUE]
            # Create warning message
            Rows <- nrow(Simb)
            if (Rows > 2000000) {
              message(paste("There are", Rows, "simulations. This may take awhile. Consider using nsim to reduce the number of simulations."))
            }
            Simb$QI <- Simb$hazard * Simb$HR
            if (!('strata' %in% names(Simb))) {
              Simb <- Simb[, list(SimID, time, QI, HRValue)]
            } else if ('strata' %in% names(Simb)) {
              Simb <- Simb[, list(SimID, time, QI, HRValue, strata)]
            }
            Simb <- data.frame(Simb)
        }
    }

    # If the user wants to calculate Hazard Rates using means for fitting all
    #covariates other than b.
    else if (isTRUE(means)) {
        if (is.null(X1) | is.null(X2)) {
            stop("For Hazard Rates, both X1 and X2 should be specified.",
                 call. = FALSE)
        }
        if (length(X1) != 1 | length(X2) != 1) {
            stop("For coxsimInteract only one value of X1 and one value of X2 can be specified.",
                call. = FALSE)
        }

        Xs <- data.frame(X1, X2)
        Xs$HRValue <- paste0(Xs$X1, ", ", Xs$X2)

        # Set all values of b at means for data used in the analysis
        NotB <- setdiff(names(DrawnDF), c(b1, b2, binter))
        MeanValues <- data.frame(obj$means)
        FittedMeans <- function(Z) {
            ID <- 1:nsim
            Temp <- data.frame(ID)
            for (i in Z) {
                BarValue <- MeanValues[i, ]
                DrawnCoef <- DrawnDF[, i]
                FittedCoef <- outer(DrawnCoef, BarValue)
                FCMolten <- data.frame(melt(FittedCoef))
                Temp <- cbind(Temp, FCMolten[,3])
            }
            Names <- c("ID", Z)
            names(Temp) <- Names
            Temp <- Temp[, -1]
            return(Temp)
        }
        FittedComb <- data.frame(FittedMeans(NotB))
        ExpandFC <- do.call(rbind, rep(list(FittedComb), nrow(Xs)))

        # Set fitted values for X1 and X2
        Simb <- data.frame(DrawnDF[, NamesInt])

        Simb <- merge(Simb, Xs)
        Simb$PreHR <- (Simb$X1 * Simb[, 1]) + (Simb$X2 * Simb[, 2]) +
                      (Simb$X1 * Simb$X2 * Simb[, 3])

        Simb <- cbind(Simb, ExpandFC)
        Simb$Sum <- rowSums(Simb[, c(7, 8)])
        Simb$HR <- exp(Simb$Sum)
        Simb <- Simb[, c("HRValue", "HR")]

        bfit <- basehaz(obj)
        bfit$FakeID <- 1
        Simb$FakeID <- 1
        bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
        SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
        Simb <- SimbDT[bfitDT, allow.cartesian = TRUE]
        # Create warning message
        Rows <- nrow(Simb)
        if (Rows > 2000000) {
          message(paste("There are", Rows,
                  "simulations. This may take awhile. Consider using nsim to reduce the number of simulations."))
        }
        Simb$QI <- Simb$hazard * Simb$HR
        if (!('strata' %in% names(Simb))) {
          Simb <- Simb[, list(time, QI, HRValue)]
        } else if ('strata' %in% names(Simb)) {
          Simb <- Simb[, list(time, QI, HRValue, strata)]
        }
        Simb <- data.frame(Simb)
    }

    # Drop simulations outside of 'confidence bounds'
    if (qi == "First Difference" | qi == "Hazard Ratio") {
        SubVar <- "X1"
    } else if (qi == "Marginal Effect") {
        SubVar <- "X2"
    } else if (qi == "Hazard Rate") {
        SubVar <- "time"
    }

    # Drop simulations outside of the middle
    SimbPerc <- IntervalConstrict(Simb = Simb, SubVar = SubVar,
                                  qi = qi, spin = spin, ci = ci,
                                  extremesDrop = extremesDrop)

    # Final clean up
    if (qi == "Hazard Rate" & !isTRUE(means)) {
        if (!('strata' %in% names(obj))) {
            SimbPercSub <- data.frame(SimbPerc$SimID,
                                      SimbPerc$time, SimbPerc$QI,
                                      SimbPerc$HRValue)
            names(SimbPercSub) <- c("SimID", "Time", "HRate",
                                    "HRValue")
        } else if ('strata' %in% names(SimbPerc)) {
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$time,
                        SimbPerc$QI, SimbPerc$strata,
                        SimbPerc$HRValue)
        names(SimbPercSub) <- c("SimID", "Time", "HRate",
                                "Strata", "HRValue")
        }
    } else if (qi == "Hazard Rate" & isTRUE(means)) {
        if (!('strata' %in% names(obj))) {
            SimbPercSub <- data.frame(SimbPerc$time, SimbPerc$QI,
                                      SimbPerc$HRValue)
            names(SimbPercSub) <- c("Time", "HRate",
                                    "HRValue")
        } else if ('strata' %in% names(SimbPerc)) {
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$time,
                        SimbPerc$QI, SimbPerc$strata,
                        SimbPerc$HRValue)
        names(SimbPercSub) <- c("SimID", "Time", "HRate",
                                "Strata", "HRValue")
        }
    } else if (qi == "Hazard Ratio") {
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$X1, SimbPerc$X2,
                                  SimbPerc$QI, SimbPerc$Comparison)
        names(SimbPercSub) <- c("SimID", "X1", "X2", "QI", "Comparison")
    } else if (qi == "Marginal Effect") {
        spalette <- NULL
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$X2, SimbPerc$QI)
        names(SimbPercSub) <- c("SimID", "X2", "QI")
    } else if (qi == "First Difference") {
        colour <- NULL
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$X1,
                                  SimbPerc$X2, SimbPerc$QI)
        names(SimbPercSub) <- c("SimID", "X1", "X2", "QI")
    }

    # Add in distribution of b
    rug_names <- names(model.frame(obj)) %in% c(b1, b2)
    rug <- model.frame(obj)[, rug_names]
    out <- list(sims = SimbPercSub, rug = rug)

    class(out) <- c("siminteract", qi, "coxsim")
    attr(out, 'b1') <- b1
    attr(out, 'b2') <- b2
    out
}
