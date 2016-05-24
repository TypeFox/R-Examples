#' Plot simulated penalised spline hazards from Cox Proportional Hazards Models
#'
#' \code{simGG.simspline} uses ggplot2 to plot
#' quantities of interest from \code{simspline} objects, including relative
#' hazards, first differences, hazard ratios, and hazard rates.
#'
#' @param obj a \code{simspline} class object
#' @param SmoothSpline logical whether or not to fit the simulations with
#' smoothing splines. Creates a smoother plot. See \code{\link{smooth.spline}}
#' for more information. Note: currently the degrees of freedom are set at 10.
#' @param FacetTime a numeric vector of points in time where you would like to
#' plot Hazard Rates in a facet grid. Only relevant if
#' \code{qi == 'Hazard Rate'}. Note: the values of Facet Time must exactly
#' match values of the \code{time} element of \code{obj}.
#' @param from numeric time to start the plot from. Only relevant if
#' \code{qi = "Hazard Rate"}.
#' @param to numeric time to plot to. Only relevant if
#' \code{qi = "Hazard Rate"}.
#' @param rug logical indicating whether or not to include a rug plot showing
#' the distribution of values in the sample used to estimate the \code{coxph}
#' model. Only relevant when the quantity of interest is not
#' \code{"Hazard Rate"}.
#' @param rug_position character string. The position adjustment to use for
#' overlapping points in the rug plot. Use \code{"jitter"} to jitter the points.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of
#' \code{qi}.
#' @param zlab a label for the plot's z-axis. Only relevant if
#' \code{qi = "Hazard Rate"} and \code{FacetTime == NULL}.
#' @param title the plot's main title.
#' @param method what type of smoothing method to use to summarize the center
#' of the simulation distribution.
#' @param lcolour character string colour of the smoothing line. The default is
#' hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if
#' \code{qi = "Relative Hazard"} or \code{qi = "First Difference"}.
#' @param lsize size of the smoothing line. Default is 1. See
#' \code{ggplot2}.
#' @param pcolour character string colour of the simulated points or ribbons
#' (when there are not multiple sets of simulations). Default is hexadecimal
#' colour \code{pcolour = '#A6CEE3'}. Only relevant if
#' \code{qi = "Relative Hazard"} or \code{qi = "First Difference"} or
#' \code{qi = "Hazard Rate"} with facets.
#' @param psize size of the plotted simulation points. Default is
#' \code{psize = 1}. See \code{ggplot2}.
#' @param alpha numeric. Alpha (e.g. transparency) for the points, lines, or
#' ribbons. Default is \code{alpha = 0.2}. See \code{ggplot2}. Note, if
#' \code{type = "lines"} or \code{type = "points"} then \code{alpah} sets the
#' maximum value per line or point at the center of the distribution. Lines or
#' points further from the center are more transparent the further they get
#' from the middle.
#' @param type character string. Specifies how to plot the simulations. Can be
#' \code{points}, \code{lines}, or \code{ribbons}. If points then each
#' simulation value will be plotted. If \code{lines} is chosen then each
#' simulation is plotted using a different line. Note: any simulation with a
#' value along its length that is outside of the specified central interval
#' will be dropped. This is to create a smooth plot. If \code{type = "ribbons"}
#' a plot will be created with shaded areas ('ribbons') for the minimum and
#' maximum simulation values (i.e. the middle interval set with \code{qi} in
#' \code{\link{coxsimSpline}}) as well as the central 50 percent of this area.
#' It also plots a line for the median value of the full area, so values in
#' \code{method} are ignored. One of the key advantages of using ribbons
#' rather than points is that it creates plots with smaller file sizes.
#' @param ... Additional arguments. (Currently ignored.)
#'
#' @return a \code{gg} \code{ggplot} class object.
#'
#' @details Uses \code{ggplot2} to plot the quantities of
#' interest from \code{simspline} objects, including relative hazards, first
#' differences, hazard ratios, and hazard rates. If currently does not support
#' hazard rates for multiple strata.
#'
#' You can to plot hazard rates for a range of values of \code{Xj} in two
#' dimensional plots at specific points in time. Each plot is arranged in a
#' facet grid.
#'
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no
#' effect, for time-varying hazard ratio graphs. No line is created for hazard
#' rates.
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' # From Keele (2010) replication data
#' M1 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 +
#'                 acutediz + hosp01  + pspline(hospdisc, df = 4) +
#'                 pspline(hhosleng, df = 4) + mandiz01 + femdiz01 +
#'                 peddiz01 + orphdum + natreg + vandavg3 + wpnoavg3 +
#'                 pspline(condavg3, df = 4) + pspline(orderent, df = 4) +
#'                 pspline(stafcder, df = 4), data = CarpenterFdaData)
#'
#' # Simulate Relative Hazards for orderent
#' Sim1 <- coxsimSpline(M1, bspline = "pspline(stafcder, df = 4)",
#'                     bdata = CarpenterFdaData$stafcder,
#'                     qi = "Hazard Ratio",
#'                     Xj = seq(1100, 1700, by = 10),
#'                     Xl = seq(1099, 1699, by = 10), spin = TRUE, nsim = 100)
#'
#' # Plot relative hazard
#' simGG(Sim1, alpha = 0.5)
#'
#' \dontrun{
#' # Simulate Hazard Rate for orderent
#' Sim2 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
#'                     bdata = CarpenterFdaData$orderent,
#'                     qi = "Hazard Rate",
#'                     Xj = seq(1, 30, by = 2), ci = 0.9, nsim = 10)
#'
#' # Create a time grid plot
#' # Find all points in time where baseline hazard was found
#' unique(Sim2$sims$Time)
#'
#' # Round time values so they can be exactly matched with FacetTime
#' Sim2$sims$Time <- round(Sim2$sims$Time, digits = 2)
#'
#' # Create plot
#' simGG(Sim2, FacetTime = c(6.21, 25.68, 100.64, 202.36),
#'        type = 'ribbons', alpha = 0.5)
#'
#' # Simulated Fitted Values of stafcder
#' Sim3 <- coxsimSpline(M1, bspline = "pspline(stafcder, df = 4)",
#'                     bdata = CarpenterFdaData$stafcder,
#'                     qi = "Hazard Ratio",
#'                     Xj = seq(1100, 1700, by = 10),
#'                     Xl = seq(1099, 1699, by = 10), ci = 0.90)
#'
#' # Plot simulated Hazard Ratios
#' simGG(Sim3, xlab = "\nFDA Drug Review Staff", type = 'lines', alpha = 0.2)
#' simGG(Sim3, xlab = "\nFDA Drug Review Staff", alpha = 0.2,
#'       SmoothSpline = TRUE, type = 'points')
#' }
#'
#' @seealso \code{\link{coxsimLinear}}, \code{\link{simGG.simtvc}},
#' \code{ggplot2}
#'
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' @import ggplot2
#' @import mgcv
#'
#' @method simGG simspline
#' @export

simGG.simspline <- function(obj, SmoothSpline = TRUE, FacetTime = NULL,
                            from = NULL, to = NULL,
                            rug = TRUE, rug_position = "identity",
                            xlab = NULL, ylab = NULL,
                            zlab = NULL, title = NULL, method = "auto",
                            lcolour = "#2B8CBE", lsize = 1, pcolour = "#A6CEE3",
                            psize = 1, alpha = 0.2, type = "ribbons", ...)
{
    Time <- Xj <- QI <- Lower50 <- Upper50 <- Min <- Max <- Median <-
    SimID <- xaxis <- NULL
    if (!inherits(obj, "simspline")){
        stop("must be a simspline object", call. = FALSE)
    }
    if (type == 'ribbons' & method != "auto"){
      message("The method argument is ignored if type = 'ribbons'. Central tendency summarised with the median.")
    }
    if (type == 'lines' & !isTRUE(SmoothSpline)){
        message(paste0('The resulting plot may look strange.',
                '\nI suggest using SmoothSpline = TRUE if type = "lines".'))
    }

    # Find quantity of interest
    qi <- class(obj)[[2]]

    if (is.null(FacetTime) & qi == 'Hazard Rate') {
        stop('FacetTime must be specified with hazard rates. scatter3d no longer
             supported.',
             .call = FALSE)
    }

    # Create y-axis label
    if (is.null(ylab)) ylab <- paste(qi, "\n")

    # Create x-axis label
    if (qi != "Hazard Rate"){
        if (is.null(xlab)) xlab <- paste("\n", attr(obj, "xaxis"))

        # Extract rug values
        rugger <- rugExtract(obj)
    }
    # Convert obj to data frame
    obj <- as.data.frame(obj)

    # Drop simulations that include outliers
    obj <- OutlierDrop(obj)

    # Smooth simulations if SmoothSpline = TRUE
    if (isTRUE(SmoothSpline)){
        if (qi != 'Hazard Rate'){
            obj <- SmoothSimulations(obj)
        } else if (qi == 'Hazard Rate'){
            obj <- SmoothSimulations(obj, xaxis = 'Time')
        }
    }

    # Alpha gradient based on percentile in the distribution
    if (type != 'ribbons' & qi != 'Hazard Rate'){
        obj <- PercRank(obj, xaxis = 'Xj')
    } else if (type != 'ribbons' & qi == 'Hazard Rate'){
        obj <- PercRank(obj, xaxis = 'Time', yaxis = 'HRate')
    }

    # Constrict time period to plot for hazard rate
    if (qi == "Hazard Rate"){
        if (is.null(xlab)) xlab <- '\nTime'
        if (!is.null(from)){
            obj <- subset(obj, Time >= from)
        }
        if (!is.null(to)){
            obj <- subset(obj, Time <= to)
        }
    }

    # Plots points
    if (type == 'points'){
        if (qi == "First Difference"){
            p <- ggplot(obj, aes(Xj, QI)) +
                    geom_point(shape = 21, aes(alpha = PercRank), size = psize,
                               colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                                color = lcolour) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
            p <- ggplot(obj, aes(Xj, QI)) +
                    geom_point(shape = 21,  aes(alpha = PercRank), size = psize,
                                                colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                                color = lcolour) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "Hazard Rate" & !is.null(FacetTime)){
            objSub <- SubsetTime(f = FacetTime, Temps = obj)
            p <- ggplot(objSub, aes(Xj, QI)) +
                    geom_point(shape = 21,  aes(alpha = PercRank), size = psize,
                                colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                                color = lcolour) +
                    facet_grid(.~Time) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE) +
                    guides(colour = guide_legend(override.aes =
                        list(alpha = 1)))
        }
    }
    # Plots lines
    else if (type == 'lines'){
        if (qi == "First Difference"){
            p <- ggplot(obj, aes(Xj, QI)) +
                    geom_line(aes(group = SimID, alpha = PercRank),
                        size = psize, colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                        color = lcolour) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
            p <- ggplot(obj, aes(Xj, QI)) +
                    geom_line(aes(group = SimID, alpha = PercRank),
                        size = psize, colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                        color = lcolour) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "Hazard Rate" & !is.null(FacetTime)){
            objSub <- SubsetTime(f = FacetTime, Temps = obj)
            p <- ggplot(objSub, aes(Xj, QI)) +
                    geom_line(aes(group = SimID, alpha = PercRank),
                        size = psize, colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                        color = lcolour) +
                    facet_grid(.~Time) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE) +
                    guides(colour = guide_legend(override.aes =
                        list(alpha = 1)))
        }
    }
    # Plots ribbons
    else if (type == 'ribbons'){
        suppressWarnings(
        if (qi == "First Difference"){
            obj <- MinMaxLines(df = obj)
            .e <- environment()
            p <- ggplot(obj, aes(Xj, Median), environment = .e) +
                    geom_line(size = lsize, colour = lcolour) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, fill = pcolour) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                        fill = pcolour) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    guides(colour = guide_legend(override.aes =
                        list(alpha = 1)))
        } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
            obj <- MinMaxLines(df = obj)
            .e <- environment()
            p <- ggplot(obj, aes(Xj, Median), environment = .e) +
                    geom_line(size = lsize, colour = lcolour) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, fill = pcolour) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                        fill = pcolour) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    guides(colour = guide_legend(override.aes =
                        list(alpha = 1)))
        } else if (qi == "Hazard Rate" & !is.null(FacetTime)){
            objSub <- SubsetTime(f = FacetTime, Temps = obj)
            objSub <- MinMaxLines(df = objSub, byVars = c("Xj", "Time"))
            .e <- environment()
            p <- ggplot(objSub, aes(Xj, Median), environment = .e) +
                    geom_line(size = lsize, colour = lcolour) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, fill = pcolour) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                        fill = pcolour) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    facet_grid(.~Time) +
                    guides(colour = guide_legend(override.aes =
                        list(alpha = 1)))
        }
        )
    }
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title) + theme_bw(base_size = 15)

    if (isTRUE(rug) & qi != 'Hazard Rate'){
        p <- p + geom_rug(data = rugger, aes(x = xaxis, y = QI), sides = "b",
                    position = rug_position, colour = pcolour)
    }
    return(p)
}
