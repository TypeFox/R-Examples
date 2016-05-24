#' Simulate quantities of interest for a range of values for a polynomial
#' nonlinear effect from Cox Proportional Hazards models
#'
#' \code{coxsimPoly} simulates quantities of interest for polynomial covariate
#' effects estimated from Cox Proportional Hazards models. These can be plotted
#' with \code{\link{simGG}}.
#' @param obj a \code{\link{coxph}} class fitted model object with a polynomial
#'  coefficient. These can be plotted with \code{\link{simGG}}.
#' @param b character string name of the coefficient you would like to simulate.
#' To find the quantity of interest using only the polynomial and not the
#' polynomial + the linear terms enter the polynomial created using
#' \code{\link{I}}, e.g. \code{I(natreg^2)} as a string.
#' @param qi quantity of interest to simulate. Values can be
#' \code{"Relative Hazard"}, \code{"First Difference"}, \code{"Hazard Ratio"},
#' and \code{"Hazard Rate"}. The default is \code{qi = "Relative Hazard"}. If
#' \code{qi = "Hazard Rate"} and the \code{coxph} model has strata, then hazard
#' rates for each strata will also be calculated.
#' @param pow numeric polynomial used in \code{coxph}.
#' @param Xj numeric vector of fitted values for \code{b} to simulate for.
#' @param Xl numeric vector of values to compare \code{Xj} to. If \code{NULL},
#' then it is authomatically set to 0.
#' @param nsim the number of simulations to run per value of \code{Xj}. Default
#' is \code{nsim = 1000}.
#' @param ci the proportion of simulations to keep. The default is
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
#' @return a \code{simpoly}, \code{coxsim} object
#'
#' @details Simulates quantities of interest for polynomial covariate effects.
#' For example if a nonlinear effect is modeled with a second order
#' polynomial--i.e. \eqn{\beta_{1}x_{i} + \beta_{2}x_{i}^{2}}{\beta[1]x[i] +
#' \beta[2]x[i]^2}--we can draw \eqn{n} simulations from the
#' multivariate normal distribution for both \eqn{\beta_{1}}{\beta[1]} and
#' \eqn{\beta_{2}}{\beta[2]}. Then we simply calculate quantities of interest
#' for a range of values and plot the results as before. For example, we find
#' the first difference for a second order polynomial with:
#' \deqn{\%\triangle h_{i}(t) = (\mathrm{e}^{\beta_{1}x_{j-1} +
#' \beta_{2}x_{j-l}^{2}} - 1) * 100}{FD(h[i](t)) = exp(\beta[1]x[j-1] +
#' \beta[2]x[j-l]^2) - 1) * 100}
#' where \eqn{x_{j-l} = x_{j} - x_{l}}{x[j-l] = x[j] - x[l]}.
#'
#' Note, you must use \code{\link{I}} to create the polynomials.
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 +
#'             acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 +
#'             peddiz01 + orphdum + natreg + I(natreg^2) +
#'             I(natreg^3) + vandavg3 + wpnoavg3 +
#'             condavg3 + orderent + stafcder, data = CarpenterFdaData)
#'
#' # Simulate simpoly First Difference
#' Sim1 <- coxsimPoly(M1, b = "natreg", qi = "First Difference",
#'                 pow = 3, Xj = seq(1, 150, by = 5), nsim = 100)
#'
#' \dontrun{
#' # Simulate simpoly Hazard Ratio with spin probibility interval
#' Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio",
#'               pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE)
#' }
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' Keele, Luke. 2010. ''Proportionally Difficult: Testing for
#' Nonproportional Hazards in Cox Models.'' Political Analysis 18(2): 189-205.
#'
#' Carpenter, Daniel P. 2002. ''Groups, the Media, Agency Waiting Costs, and
#' FDA Drug Approval.'' American Journal of Political Science 46(3): 490-505.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. ''Making the Most of
#' Statistical Analyses: Improving Interpretation and Presentation.'' American
#' Journal of Political Science 44(2): 347-61.
#'
#' Liu, Ying, Andrew Gelman, and Tian Zheng. 2013. ''Simulation-Efficient
#' Shortest Probability Intervals.'' Arvix.
#' \url{http://arxiv.org/pdf/1302.2142v1.pdf}.
#'
#' @seealso \code{\link{simGG.simpoly}}, \code{\link{survival}},
#' \code{\link{strata}}, and \code{\link{coxph}}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats vcov model.frame
#' @importFrom survival basehaz
#' @importFrom dplyr inner_join as_data_frame
#' @export

coxsimPoly <- function(obj, b = NULL, qi = "Relative Hazard", pow = 2,
                       Xj = NULL, Xl = NULL, nsim = 1000, ci = 0.95,
                       spin = FALSE, extremesDrop = TRUE)
{
    strata <- QI <- SimID <- time <- NULL
    # Ensure that qi is valid
    qiOpts <- c("Relative Hazard", "First Difference", "Hazard Rate",
             "Hazard Ratio")
    TestqiOpts <- qi %in% qiOpts
    if (!isTRUE(TestqiOpts)){
        stop("Invalid qi type.\nqi must be 'Relative Hazard', 'Hazard Rate', 'First Difference', or 'Hazard Ratio'.",
        call. = FALSE)
    }
    # Ensure that b is declared
    if (is.null(b)) stop('Need to declare b.', call. = FALSE)

    # Find base variable if polynomial entered
    PolyOnly <- grepl(pattern = 'I\\(.*\\^', x = b)
    if (isTRUE(PolyOnly)){
        message('Simulations of the quantity of interest will not include the linear term, if any exist in the model.\n')
        b <- gsub(pattern = '^I\\(', replacement = '', x = b)
        b <- gsub(pattern = '\\^.*$', replacement = '', x = b)
    }

    # Find X_{jl}
    if (length(Xj) != length(Xl) & !is.null(Xl)){
        stop("Xj and Xl must be the same length.", call. = FALSE)
    } else if (is.null(Xl)) {
        message("All Xl set at 0.")
        Xjl <- Xj
    } else {
        Xbound <- cbind(Xj, Xl)
        Xjl <- Xbound[, 1] - Xbound[, 2]
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

    # Subset data frame to only include polynomial constitutive terms.
    bpos <- match(b, dfn)
    NamesLoc <- function(p){
        Temp <- paste0("I.", b, ".", p, ".")
        match(Temp, dfn)
    }
    pows <- as.numeric(2:pow)
    NamePow <- sapply(pows, NamesLoc, simplify = TRUE)

    if (!isTRUE(PolyOnly)) NamePow <- c(bpos, NamePow)

    Drawn <- data.frame(Drawn[, NamePow])
    VNames <- names(Drawn)
    powFull <- as.numeric(1:pow)

    # Function to Multiply covariates by polynomials
    Fitted <- function(VN, x, p){
        Temp <- outer(Drawn[, VN], x^p)
        TempDF <- data.frame(melt(Temp))
        TempDF <- TempDF[, 'value']
        TempDF
    }

    Simb <- data.frame()
    for (i in Xjl){
        TempComb <- mapply(Fitted, VN = VNames, x = i, p = powFull)
        TempComb <- data.frame(TempComb)
        TempComb$Xjl <- i
        TempComb$SimID <- SimID
        Simb <- rbind(Simb, TempComb)
    }

    # Create combined quantities of interest
    if (length(NamePow) > 1){
        UnExp <- rowSums(Simb[, VNames])
    }
    else if (length(NamePow) == 1) UnExp = Simb[, VNames]

    if (qi == "Relative Hazard"){
        Simb$QI <- exp(UnExp)
    } else if (qi == "Hazard Ratio"){
        Simb$QI <- exp(UnExp)
    }
    else if (qi == "First Difference"){
        Simb$QI <- (exp(UnExp) - 1) * 100
    }
    else if (qi == "Hazard Rate"){
        Simb$HR <- exp(UnExp)
        bfit <- basehaz(obj)
        bfit$FakeID <- 1
        Simb$FakeID <- 1
        bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
        SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
        Simb <- SimbDT[bfitDT, allow.cartesian = TRUE]
        # Create warning message
        Rows <- nrow(Simb)
        if (Rows > 2000000){
            message(paste("There are", Rows,
            "simulations. This may take awhile. Consider using nsim to reduce the number of simulations."))
        }
        Simb$QI <- Simb$hazard * Simb$HR
        if (!('strata' %in% names(Simb))){
            Simb <- Simb[, list(SimID, time, Xjl, QI)]
        } else if ('strata' %in% names(Simb)){
        Simb <- Simb[, list(SimID, time, Xjl, QI, strata)]
        }
        Simb <- data.frame(Simb)
    }

    # Drop simulations outside of 'confidence bounds'
    if (qi != "Hazard Rate"){
        SubVar <- "Xjl"
    } else if (qi == "Hazard Rate"){
        names(Simb)[names(Simb) == "Xjl"] <- "HRValue"
        SubVar <- c("SimID", "time", "HRValue")
    }

    # Drop simulations outside of the middle
    SimbPerc <- IntervalConstrict(Simb = Simb, SubVar = SubVar,
                                qi = qi, spin = spin, ci = ci,
                                extremesDrop = extremesDrop)

    # Clean up
    if (qi == "Hazard Rate"){
        if (!('strata' %in% names(SimbPerc))){
            SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$time,
                                    SimbPerc$QI, SimbPerc$HRValue)
            names(SimbPercSub) <- c("SimID", "Time", "HRate", "HRValue")
        } else if ('strata' %in% names(SimbPerc)) {
            SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$time,
                                SimbPerc$QI, SimbPerc$strata, SimbPerc$HRValue)
            names(SimbPercSub) <- c("SimID", "Time", "HRate", "Strata",
                                    "HRValue")
        }
    } else if (qi == "Hazard Ratio" | qi == "Relative Hazard" |
            qi == "First Difference"){
        merger_xj <- data.frame(Xj = Xj, Xjl = Xjl)
        SimbPerc <- inner_join(SimbPerc, merger_xj, by = 'Xjl') %>%
                        as_data_frame
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$Xj, SimbPerc$QI)
        names(SimbPercSub) <- c("SimID", "Xj", "QI")
    }

    # Add in distribution of b
    rug <- model.frame(obj)[, b]
    out <- list(sims = SimbPercSub, rug = rug)

    class(out) <- c("simpoly", qi, "coxsim")
    attr(out, 'xaxis') <- b
    out
}
