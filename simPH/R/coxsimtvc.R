#' Simulate time-interactive quantities of interest from Cox Proportional
#' Hazards models
#'
#' \code{coxsimtvc} simulates time-interactive relative hazards, first
#' differences, and hazard ratios from models estimated with \code{\link{coxph}}
#' using the multivariate normal distribution. These can be plotted with
#' \code{\link{simGG}}.
#' @param obj a \code{\link{coxph}} fitted model object with a time interaction.
#' @param b the non-time interacted variable's name.
#' @param btvc the time interacted variable's name.
#' @param qi character string indicating what quantity of interest you would
#' like to calculate. Can be \code{'Relative Hazard'},
#' \code{'First Difference'}, \code{'Hazard Ratio'}, \code{'Hazard Rate'}.
#' Default is \code{qi = 'Relative Hazard'}. If \code{qi = 'First Difference'}
#' or \code{qi = 'Hazard Ratio'} then you can set \code{Xj} and \code{Xl}.
#' @param Xj numeric vector of fitted values for \code{b}. Must be the same
#' length as \code{Xl} or \code{Xl} must be \code{NULL}.
#' @param Xl numeric vector of fitted values for Xl. Must be the same length as
#' \code{Xj}. Only applies if \code{qi = 'First Difference'} or
#' \code{qi = 'Hazard Ratio'}.
#' @param nsim the number of simulations to run per point in time. Default is
#' \code{nsim = 1000}.
#' @param tfun function of time that btvc was multiplied by. Default is
#' "linear". It can also be "log" (natural log) and "power". If
#' \code{tfun = "power"} then the pow argument needs to be specified also.
#' @param pow if \code{tfun = "power"}, then use pow to specify what power the
#' time interaction was raised to.
#' @param from point in time from when to begin simulating coefficient values
#' @param to point in time to stop simulating coefficient values.
#' @param by time intervals by which to simulate coefficient values.
#' @param ci the proportion of simulations to keep. The default is
#' \code{ci = 0.95}, i.e. keep the middle 95 percent. If \code{spin = TRUE}
#' then \code{ci} is the confidence level of the shortest probability interval.
#' Any value from 0 through 1 may be used.
#' @param spin logical, whether or not to keep only the shortest probability
#' interval rather than the middle simulations. Currently not supported for
#' hazard rates.
#' @param extremesDrop logical whether or not to drop simulated quantity of
#' interest values that are \code{Inf}, \code{NA}, \code{NaN} and
#' \eqn{> 1000000} for \code{spin = FALSE} or \eqn{> 800} for
#' \code{spin = TRUE}.
#' These values are difficult to plot \code{\link{simGG}} and may prevent
#' \code{spin} from finding the central interval.
#'
#' @return a \code{simtvc} object
#'
#' @details Simulates time-varying relative hazards, first differences, and
#' hazard ratios using parameter estimates from \code{coxph} models. Can also
#' simulate hazard rates for multiple strata.
#'
#' Relative hazards are found using:
#' \deqn{RH = e^{\beta_{1} + \beta_{2}f(t)}}{RH = exp(\beta[1] + \beta[2] f(t))}
#' where \eqn{f(t)} is the function of time.
#'
#' First differences are found using:
#' \deqn{FD = (e^{(X_{j} - X_{l}) (\beta_{1} + \beta_{2}f(t))} - 1) * 100}{FD =
#' (exp((X[j] - X[l])(\beta[1] + \beta[2]f(t)) - 1) * 100}
#' where \eqn{X_{j}}{X[j]} and \eqn{X_{l}}{X[l]} are some values of \eqn{X} to
#' contrast.
#'
#' Hazard ratios are calculated using:
#' \deqn{FD = e^{(X_{j} - X_{l}) (\beta_{1} + \beta_{2}f(t))}}{FD = exp((X[j] - X[l])(\beta[1] + \beta[2]f(t))}
#' When simulating non-stratifed time-varying harzards \code{coxsimtvc} uses
#' the point estimates for a given coefficient \eqn{\hat{\beta}_{x}}{hat \beta[x]} and its time interaction \eqn{\hat{\beta}_{xt}}{hat \beta[xt]}
#' along with the variance matrix (\eqn{\hat{V}(\hat{\beta})}{hat V(hat \beta)})
#' estimated from a \code{coxph} model. These are used to draw values of
#' \eqn{\beta_{1}}{\beta[1]} and \eqn{\beta_{2}}{\beta[2]} from the
#' multivariate normal distribution \eqn{N(\hat{\beta},\: \hat{V}(\hat{\beta}))}{N(hat \beta, hat V (hat \beta))}.
#'
#' When simulating stratified time-varying hazard rates \eqn{H} for a given
#' strata \eqn{k}, \code{coxsimtvc} uses:
#' \deqn{H_{kxt} = \hat{\beta}_{k0t}e^{\hat{\beta_{1}} + \beta_{2}f(t)}}{H_{kxt} = hat \beta[k0t]exp(hat \beta[1] + \beta[2]f(t))}
#' The resulting simulation values can be plotted using \code{\link{simGG}}.
#'
#' @examples
#' \dontrun{
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Expand data (not run to speed processing time, but should be run)
#' GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
#'                      Time = 'begin', Time2 = 'end', event = 'event')
#'
#' # Create time interactions
#' BaseVars <- c('qmv', 'backlog', 'coop', 'codec', 'qmvpostsea', 'thatcher')
#' GolubEUPData <- tvc(GolubEUPData, b = BaseVars, tvar = 'end', tfun = 'log')
#'
#' # Run Cox PH Model
#' M1 <- coxph(Surv(begin, end, event) ~ qmv + qmvpostsea + qmvpostteu +
#'                 coop + codec + eu9 + eu10 + eu12 + eu15 + thatcher +
#'                 agenda + backlog + qmv_log + qmvpostsea_log + coop_log +
#'                 codec_log + thatcher_log + backlog_log,
#'             data = GolubEUPData, ties = "efron")
#'
#' # Create simtvc object for Relative Hazard
#' Sim1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "qmv_log",
#'                    tfun = "log", from = 80, to = 2000,
#'                    Xj = 1, by = 15, ci = 0.99, nsim = 100)
#'
#' # Create simtvc object for First Difference
#' Sim2 <- coxsimtvc(obj = M1, b = "qmv", btvc = "qmv_log",
#'                  qi = "First Difference", Xj = 1,
#'                  tfun = "log", from = 80, to = 2000,
#'                  by = 15, ci = 0.95)
#'
#' # Create simtvc object for Hazard Ratio
#' Sim3 <- coxsimtvc(obj = M1, b = "backlog", btvc = "backlog_log",
#'                   qi = "Hazard Ratio", Xj = c(191, 229),
#'                   Xl = c(0, 0),
#'                   tfun = "log", from = 80, to = 2000,
#'                   by = 15, ci = 0.5)
#' }
#'
#' @seealso \code{\link{simGG}}, \code{\link{survival}}, \code{\link{strata}},
#' and \code{\link{coxph}}
#'
#' @import data.table
#' @importFrom stats vcov
#' @importFrom survival basehaz
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' Golub, Jonathan, and Bernard Steunenberg. 2007. ''How Time
#' Affects EU Decision-Making.'' European Union Politics 8(4): 555-66.
#'
#' Licht, Amanda A. 2011. ''Change Comes with Time: Substantive Interpretation
#' of Nonproportional Hazards in Event History Analysis.'' Political Analysis
#' 19: 227-43.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. ''Making the Most of
#' Statistical Analyses: Improving Interpretation and Presentation.'' American
#' Journal of Political Science 44(2): 347-61.
#'
#' Liu, Ying, Andrew Gelman, and Tian Zheng. 2013. ''Simulation-Efficient
#' Shortest Probability Intervals.'' Arvix.
#' \url{http://arxiv.org/pdf/1302.2142v1.pdf}.

coxsimtvc <- function(obj, b, btvc, qi = "Relative Hazard", Xj = NULL,
                      Xl = NULL, tfun = "linear", pow = NULL, nsim = 1000,
                      from, to, by = 1, ci = 0.95, spin = FALSE,
                      extremesDrop = TRUE)
{
    HRValue <- strata <- QI <- NULL
    # Ensure that qi is valid
    qiOpts <- c("Relative Hazard", "First Difference", "Hazard Rate",
                "Hazard Ratio")
    TestqiOpts <- qi %in% qiOpts
    if (!isTRUE(TestqiOpts)){
    stop("Invalid qi type. qi must be 'Relative Hazard', 'First Difference',
        'Hazard Rate', or 'Hazard Ratio'.",
        call. = FALSE)
    }

    if (is.null(Xl) & qi != "Hazard Rate"){
        Xl <- rep(0, length(Xj))
        message("All Xl set to 0.")
    } else if (!is.null(Xl) & qi == "Relative Hazard") {
        message("All Xl set to 0.")
    }

    # Create time function
    tfunOpts <- c("linear", "log", "power")
    TestforTOpts <- tfun %in% tfunOpts
    if (!isTRUE(TestforTOpts)){
        stop("Must specify tfun as 'linear', 'log', or 'power'", call. = FALSE)
    }

    if (tfun == "linear"){
        tf <- seq(from = from, to = to, by = by)
    } else if (tfun == "log"){
        tf <- log(seq(from = from, to = to, by = by))
    } else if (tfun == "power"){
        tf <- (seq(from = from, to = to, by = by))^pow
    }

    # Create simulation ID variable
    SimID <- 1:nsim

    # Parameter estimates & Varance/Covariance matrix
    Coef <- matrix(obj$coefficients)
    VC <- vcov(obj)

    # Draw values from the multivariate normal distribution
    Drawn <- mvrnorm(n = nsim, mu = Coef, Sigma = VC)
    DrawnDF <- data.frame(Drawn)
    dfn <- names(DrawnDF)

    # Extract simulations for variables of interest
    bpos <- match(b, dfn)
    btvcpos <- match(btvc, dfn)

    Drawn <- data.frame(Drawn[, c(bpos, btvcpos)])
    Drawn$SimID <- SimID

    # Multiply time function with btvc
    TVSim <- outer(Drawn[,2], tf)
    TVSim <- data.frame(melt(TVSim))
    names(TVSim) <- c("SimID", "time", "TVC")
    time <- 1:length(tf)
    TempDF <- data.frame(time, tf)
    TVSim <- merge(TVSim, TempDF)

    # Combine with non TVC version of the variable
    TVSim <- merge(Drawn, TVSim, by = "SimID")
    TVSim$CombCoef <- TVSim[[2]] + TVSim$TVC

    # Find quantity of interest
  if (qi == "Relative Hazard"){
        Xs <- data.frame(Xj)
        names(Xs) <- c("Xj")
        Xs$Comparison <- paste(Xs[, 1])
        Simb <- merge(TVSim, Xs)
        Simb$QI <- exp(Simb$CombCoef * Simb$Xj)
    } else if (qi == "First Difference"){
        if (length(Xj) != length(Xl)){
            stop("Xj and Xl must be the same length.", call. = FALSE)
        } else {
            TVSim$QI <- exp(TVSim$CombCoef)
            Xs <- data.frame(Xj, Xl)
            Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
            Simb <- merge(TVSim, Xs)
            Simb$QI <- (exp((Simb$Xj - Simb$Xl) * Simb$CombCoef) - 1) * 100
        }
    } else if (qi == "Hazard Ratio"){
        if (length(Xj) != length(Xl)){
            stop("Xj and Xl must be the same length.", call. = FALSE)
        } else {
            Xs <- data.frame(Xj, Xl)
            Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
            Simb <- merge(TVSim, Xs)
            Simb$QI <- exp((Simb$Xj - Simb$Xl) * Simb$CombCoef)
        }
    } else if (qi == "Hazard Rate"){
        if (!is.null(Xl)) {
            Xl <- NULL
            message("Xl is ignored.")
        }
        Xs <- data.frame(Xj)
        Xs$HRValue <- paste(Xs[, 1])
        Simb <- merge(TVSim, Xs)
        Simb$HR <- exp(Simb$Xj * Simb$CombCoef)
        bfit <- basehaz(obj)
        bfit$FakeID <- 1
        Simb$FakeID <- 1
        bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
        SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
        Simb <- SimbDT[bfitDT, allow.cartesian = TRUE]
        # Create warning message
        Rows <- nrow(Simb)
        if (Rows > 2000000){
            message(paste("There are", Rows, "simulations. This may take awhile. Consider using nsim to reduce the number of simulations."))
        }
        Simb$QI <- Simb$hazard * Simb$HR
        if (!('strata' %in% names(Simb))){
            Simb <- Simb[, list(SimID, time, tf, Xj, QI, HRValue)]
        } else if ('strata' %in% names(Simb)){
            Simb <- Simb[, list(SimID, time, tf, Xj, QI, HRValue, strata)]
        }
        Simb <- data.frame(Simb)
    }

  # Drop simulations outside of 'confidence bounds'
  SubVar <- c("time", "Xj")

  # Drop simulations outside of the middle
  SimbPerc <- IntervalConstrict(Simb = Simb, SubVar = SubVar,
                                qi = qi, spin = spin, ci = ci,
                                extremesDrop = extremesDrop)

    # Create real time variable
    if (tfun == "linear"){
        SimbPerc$RealTime <- SimbPerc$tf
    } else if (tfun == "log"){
        SimbPerc$RealTime <- exp(SimbPerc$tf)
    } else if (tfun == "power"){
        SimbPerc$RealTime <- SimbPerc$tf^(1/pow)
    }

    # Final clean up
    # Subset simtvc object & create data frame of important variables
    if (qi == "Hazard Rate"){
        if (!('strata' %in% names(SimbPerc))){
            SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$RealTime,
                                SimbPerc$QI,
            SimbPerc$HRValue)
            names(SimbPercSub) <- c("SimID", "Time", "HRate", "HRValue")
        } else if ('strata' %in% names(SimbPerc)) {
            SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$RealTime,
                            SimbPerc$QI,
            SimbPerc$strata, SimbPerc$HRValue)
            names(SimbPercSub) <- c("SimID", "Time", "HRate", "Strata",
                                    "HRValue")
        }
    } else if (qi == "Hazard Ratio"){
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$RealTime, SimbPerc$QI,
                                SimbPerc$Comparison)
        names(SimbPercSub) <- c("SimID", "Time", "QI", "Comparison")
    } else if (qi == "Relative Hazard"){
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$RealTime,
            SimbPerc$QI, SimbPerc$Xj)
        names(SimbPercSub) <- c("SimID", "Time", "QI", "Xj")
    } else if (qi == "First Difference"){
        SimbPercSub <- data.frame(SimbPerc$SimID, SimbPerc$RealTime,
                            SimbPerc$QI, SimbPerc$Comparison)
        names(SimbPercSub) <- c("SimID", "Time", "QI", "Comparison")
    }
    class(SimbPercSub) <- c("simtvc", qi, "data.frame")
    SimbPercSub
}
