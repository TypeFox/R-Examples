#' Plot simulated time-interactive hazard ratios or stratified time-interactive
#' hazard rates from Cox Proportional Hazards Models
#'
#' \code{simGG.simtvc} uses ggplot2 to plot the simulated hazards from a
#' \code{simtvc} class object created by \code{\link{coxsimtvc}} using
#' ggplot2.
#' @param obj a \code{simtvc} class object
#' @param from numeric time to start the plot from.
#' @param to numeric time to plot to.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of
#' \code{qi}.
#' @param title the plot's main title.
#' @param method what type of smoothing method to use to summarize the center
#' of the simulation distribution.
#' @param spalette colour palette for when there are multiple sets of
#' comparisons to plot. Default palette is \code{"Set1"}. See
#' \code{\link{scale_colour_brewer}}.
#' @param legend specifies what type of legend to include (if applicable).
#' The default is \code{legend = "legend"}. To hide the legend use
#' \code{legend = FALSE}. See the \code{\link{discrete_scale}} for more details.
#' @param leg.name name of the legend (if applicable).
#' @param lsize size of the smoothing line. Default is 1. See
#' \code{ggplot2}.
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
#' @return a \code{gg} \code{ggplot} class object
#' @details Plots either a time-interactive hazard ratios, first differences,
#' and relative hazards, or the hazard rates for multiple strata. Currently the
#' strata legend labels need to be changed manually (see \code{revalue}
#' in the plyr package) in the \code{simtvc} object with the
#' \code{strata} component. Also, currently the x-axis tick marks and break
#' labels must be adjusted manually for non-linear functions of time.
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no
#' effect, for time-varying hazard ratio graphs. No line is created for hazard
#' rates.
#'
#' @examples
#' \dontrun{
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Expand data
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
#' # Create plot
#' simGG(Sim1, legend = FALSE)
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
#'                   tfun = "log", from = 100, to = 2000,
#'                   by = 15, ci = 0.99)
#'
#' # Create plots
#' simGG(Sim2, type = 'points')
#' simGG(Sim3, leg.name = "Comparision", from = 1200, type = 'lines')
#' }
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' Licht, Amanda A. 2011. ''Change Comes with Time: Substantive
#' Interpretation of Nonproportional Hazards in Event History Analysis.''
#' Political Analysis 19: 227-43.
#' @import ggplot2
#' @import mgcv
#'
#' @method simGG simtvc
#' @export

simGG.simtvc <- function(obj, from = NULL, to = NULL, xlab = NULL, ylab = NULL,
                        title = NULL, method = "auto", spalette = "Set1",
                        legend = "legend", leg.name = "", lsize = 1, psize = 1,
                        alpha = 0.2, type = "ribbons", ...)
{
    Time <- HRate <- HRValue <- QI <- Comparison <- Xj <- Lower50 <- Upper50 <-
    Min <- Max <- Median <- SimID <- NULL
    if (!inherits(obj, "simtvc")){
        stop("must be a simtvc object", call. = FALSE)
    }
    if (type == "ribbons" & method != "auto"){
        message("The method argument is ignored if type = 'ribbons'. Central tendency summarised with the median.")
    }

    # Find quantity of interest
    qi <- class(obj)[[2]]

    # Create y-axis label
    if (is.null(ylab)) ylab <- paste(qi, "\n")

    # Create x-axis label
    if (is.null(xlab)) xlab <- '\n Time'

    # Convert obj to data frame
    obj <- as.data.frame(obj)

    # Drop simulations that include outliers
    if (type == 'lines'){
        obj <- OutlierDrop(obj)
    }

    # Alpha gradient based on percentile in the distribution
    if (type != 'ribbons' & qi != 'Hazard Rate'){
        obj <- PercRank(obj, xaxis = 'Time')
    } else if (type != 'ribbons' & qi == 'Hazard Rate'){
        obj <- PercRank(obj, xaxis = 'Time', yaxis = 'HRate')
    }

    # Constrict time period
    if (!is.null(from)){
        obj <- subset(obj, Time >= from)
    }
    if (!is.null(to)){
        obj <- subset(obj, Time <= to)
    }

    # Plot points
    if (type == 'points'){
        if (qi == "Hazard Rate"){
            if ('strata' %in% names(obj)) {
                p <- ggplot(obj, aes(x = Time, y = HRate,
                        colour = factor(HRValue))) +
                        geom_point(aes(alpha = PercRank), size = psize) +
                        geom_smooth(method = method, size = lsize, se = FALSE) +
                        facet_grid(.~ Strata) +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                            guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            } else if (!('strata' %in% names(obj))){
                p <- ggplot(obj, aes(Time, HRate, colour = factor(HRValue))) +
                        geom_point(shape = 21, aes(alpha = PercRank),
                            size = psize) +
                        geom_smooth(method = method, size = lsize, se = FALSE) +
                        scale_colour_brewer(palette = spalette,
                                name = leg.name, guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        } else if (qi == "Hazard Ratio"){
            p <- ggplot(obj, aes(x = Time, y = QI,
                    colour = factor(Comparison))) +
                    geom_point(aes(alpha = PercRank), size = psize) +
                    geom_smooth(method = method, size = lsize, se = FALSE) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                                      guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "Relative Hazard"){
            p<- ggplot(obj, aes(x = Time, y = QI, colour = factor(Xj))) +
                    geom_point(aes(alpha = PercRank), size = psize) +
                    geom_smooth(method = method, size = lsize, se = FALSE) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                                      guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "First Difference"){
            p <- ggplot(obj, aes(Time, QI, colour = factor(Comparison))) +
                    geom_point(shape = 21, aes(alpha = PercRank), size = psize) +
                    geom_smooth(method = method, size = lsize, se = FALSE) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                                      guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
            }
    }
    # Plot lines
    else if (type == 'lines'){
        if (qi == "Hazard Rate"){
            if ('strata' %in% names(obj)) {
                p <- ggplot(obj, aes(x = Time, y = HRate,
                        colour = factor(HRValue))) +
                        geom_line(aes(group = interaction(SimID,
                            factor(HRValue)), alpha = PercRank), size = psize) +
                        geom_smooth(aes(colour = factor(HRValue)),
                            method = method, size = lsize, se = FALSE) +
                        facet_grid(.~ Strata) +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                                          guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            } else if (!('strata' %in% names(obj))){
                p <- ggplot(obj, aes(Time, HRate, colour = factor(HRValue))) +
                    geom_line(aes(group = interaction(SimID, factor(HRValue)),
                                alpha = PercRank), shape = 21, size = psize) +
                    geom_smooth(aes(colour = factor(HRValue)), method = method,
                                  size = lsize, se = FALSE) +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                                      guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
            }
        } else if (qi == "Hazard Ratio"){
            p <- ggplot(obj, aes(x = Time, y = QI, colour = factor(Comparison))) +
                    geom_line(aes(group = interaction(SimID,
                        factor(Comparison)), alpha = PercRank), size = psize) +
                    geom_smooth(aes(group = factor(Comparison)),
                        method = method, size = lsize, se = FALSE) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                        guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "Relative Hazard"){
            p <- ggplot(obj, aes(x = Time, y = QI, colour = factor(Xj))) +
                    geom_line(aes(group = interaction(SimID, factor(Xj)),
                        alpha = PercRank), size = psize) +
                    geom_smooth(aes(group = factor(Xj)), method = method,
                        size = lsize, se = FALSE) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                        guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "First Difference"){
            p <- ggplot(obj, aes(Time, QI, colour = factor(Comparison))) +
                    geom_line(aes(group = interaction(SimID,
                        factor(Comparison)), alpha = PercRank), size = psize) +
                    geom_smooth(aes(group = factor(Comparison)),
                        method = method, size = lsize, se = FALSE) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                        guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        }
    }
    # Plot ribbons
    else if (type == 'ribbons'){
        suppressWarnings(
        if (qi == "Hazard Rate"){
            if ('strata' %in% names(obj)) {
                obj <- MinMaxLines(df = obj, hr = TRUE, strata = TRUE)
                p <- ggplot(obj, aes(x = Time, y = HRate,
                        colour = factor(HRValue), fill = factor(HRValue))) +
                        geom_line(size = lsize) +
                        geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                            alpha = alpha, linetype = 0) +
                        geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                            linetype = 0) +
                        facet_grid(. ~ Strata) +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                            guide = legend) +
                        scale_fill_brewer(palette = spalette, name = leg.name,
                            guide = legend)
            } else if (!('strata' %in% names(obj))){
                obj <- MinMaxLines(df = obj, hr = TRUE)
                p <- ggplot(obj, aes(Time, Median, colour = factor(HRValue),
                        fill = factor(HRValue))) +
                        geom_line(size = lsize) +
                        geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                            alpha = alpha, linetype = 0) +
                        geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                            linetype = 0) +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                            guide = legend) +
                        scale_fill_brewer(palette = spalette, name = leg.name,
                            guide = legend)
            }
        } else if (qi == "Hazard Ratio"){
            obj <- MinMaxLines(df = obj, byVars = c("Time", "Comparison"))
            p <- ggplot(obj, aes(x = Time, y = Median,
                    colour = factor(Comparison), fill = factor(Comparison))) +
                    geom_line(size = lsize) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, linetype = 0) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                        linetype = 0) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                        guide = legend) +
                    scale_fill_brewer(palette = spalette, name = leg.name,
                        guide = legend)
        } else if (qi == "Relative Hazard"){
            obj <- MinMaxLines(df = obj, byVars = c("Time", "Xj"))
            p <- ggplot(obj, aes(x = Time, y = Median, colour = factor(Xj),
                    fill = factor(Xj))) +
                    geom_line(size = lsize) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, linetype = 0) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                        linetype = 0) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                        guide = legend) +
                    scale_fill_brewer(palette = spalette, name = leg.name,
                        guide = legend)
        } else if (qi == "First Difference"){
            obj <- MinMaxLines(df = obj, byVars = c("Time", "Comparison"))
            p <- ggplot(obj, aes(x = Time, y = Median,
                    colour = factor(Comparison), fill = factor(Comparison))) +
                    geom_line(size = lsize) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, linetype = 0) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                        linetype = 0) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                        guide = legend) +
                    scale_fill_brewer(palette = spalette, name = leg.name,
                        guide = legend)
        }
        )
    }
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title) + theme_bw(base_size = 15)

    return(p)
}
