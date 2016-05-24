#' Plot simulated polynomial quantities of interest from Cox Proportional
#' Hazards Models
#'
#' \code{simGG.simpoly} uses ggplot2 to plot simulated relative
#' quantities of interest from a \code{simpoly} class object.
#' @param obj a \code{simpoly} class object.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of
#' \code{qi}.
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
#' @param title the plot's main title.
#' @param method what type of smoothing method to use to summarize the center
#' of the simulation distribution.
#' @param spalette colour palette for when there are multiple sets of
#' comparisons to plot. Default palette is \code{"Set1"}. See
#' \code{\link{scale_colour_brewer}}.
#' @param legend specifies what type of legend to include (if applicable). The
#' default is \code{legend = "legend"}. To hide the legend use
#' \code{legend = FALSE}. See the \code{\link{discrete_scale}} for more details.
#' @param leg.name name of the legend (if applicable).
#' @param lcolour character string colour of the smoothing line. The default is
#' hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if
#' \code{qi = "First Difference"}.
#' @param lsize size of the smoothing line. Default is 1. See
#' \code{ggplot2}.
#' @param pcolour character string colour of the simulated points or ribbons
#' (when there are not multiple sets of simulations). Default is hexadecimal
#' colour \code{pcolour = '#A6CEE3'}.
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
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
#'        deathrt1 + acutediz + hosp01  + hhosleng + mandiz01 +
#'        femdiz01 + peddiz01 + orphdum + natreg +
#'        I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 +
#'        condavg3 + orderent + stafcder, data = CarpenterFdaData)
#'
#' # Simulate simpoly First Difference
#' Sim1 <- coxsimPoly(M1, b = "natreg", qi = "First Difference",
#'            pow = 3, Xj = seq(1, 150, by = 5), nsim = 100)
#'
#' # Plot simulations
#' simGG(Sim1, rug_position = 'jitter')
#'
#' \dontrun{
#' # Simulate simpoly Hazard Ratio with spin probibility interval
#' Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio",
#'           pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE,
#'           nsim = 100)
#'
#' # Plot simulations
#' simGG(Sim2, type = 'ribbons', rug_position = 'jitter')
#'
#' Sim3 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Rate",
#'            pow = 3, Xj = c(1, 150), nsim = 100)
#'
#' # Plot simulations
#' simGG(Sim3, type = 'lines')
#' }
#'
#' @details Uses ggplot2 to plot the quantities of interest from
#' \code{simpoly} objects.
#'
#' @seealso \code{\link{coxsimPoly}} and \code{ggplot2}
#'
#' @return a \code{gg} \code{ggplot} class object
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' @import ggplot2
#' @import mgcv
#'
#' @method simGG simpoly
#' @export

simGG.simpoly <- function(obj, from = NULL, to = NULL,
                        rug = TRUE, rug_position = "identity",
                        xlab = NULL, ylab = NULL,
                        title = NULL, method = "auto", spalette = "Set1",
                        legend = "legend", leg.name = "", lcolour = "#2B8CBE",
                        lsize = 1, pcolour = "#A6CEE3", psize = 1,
                        alpha = 0.2, type = "ribbons", ...)
{
    Time <- HRValue <- HRate <- Xj <- QI <- Lower50 <- Upper50 <- Min <- Max <-
    Median <- SimID <- xaxis <- NULL
    if (!inherits(obj, "simpoly")){
        stop("must be a simpoly object", call. = FALSE)
    }
    if (type == 'ribbons' & method != "auto"){
        message("The method argument is ignored if type = 'ribbons'. Central tendency summarised with the median.")
    }
    # Find quantity of interest
    qi <- class(obj)[[2]]

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
    if (type == 'lines'){
        obj <- OutlierDrop(obj)
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
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
            } else if (!('strata' %in% names(obj))){
                p <- ggplot(obj, aes(Time, HRate, colour = factor(HRValue))) +
                    geom_point(shape = 21, aes(alpha = PercRank), size = psize) +
                    geom_smooth(method = method, size = lsize, se = FALSE) +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                                        guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
            }
        } else if (qi == "First Difference"){
            p <- ggplot(obj, aes(Xj, QI)) +
                    geom_point(shape = 21, aes(alpha = PercRank),
                                size = psize, colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                              color = lcolour) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
            p <- ggplot(obj, aes(Xj, QI)) +
                    geom_point(shape = 21, aes(alpha = PercRank),
                                size = psize, colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                              color = lcolour) +
                    scale_alpha_continuous(range = c(0, alpha),
                                            guide = FALSE) +
                    geom_hline(aes(yintercept = 1), linetype = "dotted")
        }
    }
    # Plot lines
    else if (type == 'lines'){
        if (qi == "Hazard Rate"){
            if ('strata' %in% names(obj)) {
                p <- ggplot(obj, aes(x = Time, y = HRate,
                        colour = factor(HRValue))) +
                        geom_line(aes(group = interaction(SimID,
                            factor(HRValue)), alpha = PercRank),
                            size = psize) +
                        geom_smooth(aes(colour = factor(HRValue)),
                            method = method, size = lsize, se = FALSE) +
                        facet_grid(.~ Strata) +
                        scale_colour_brewer(palette = spalette,
                                name = leg.name, guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                                guide = FALSE)
            } else if (!('strata' %in% names(obj))){
                p <- ggplot(obj, aes(Time, HRate,
                        colour = factor(HRValue))) +
                        geom_line(aes(group = interaction(SimID,
                            factor(HRValue)), alpha = PercRank),
                            size = psize) +
                        geom_smooth(aes(colour = factor(HRValue)),
                            method = method, size = lsize, se = FALSE) +
                        scale_colour_brewer(palette = spalette,
                            name = leg.name, guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        } else if (qi == "First Difference"){
            p <- ggplot(obj, aes(Xj, QI)) +
                    geom_line(aes(group = SimID, alpha = PercRank),
                            size = psize, colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                              color = lcolour) +
                    geom_hline(aes(yintercept = 0), linetype = "dotted") +
                    scale_alpha_continuous(range = c(0, alpha),
                        guide = FALSE)
        } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
            p <- ggplot(obj, aes(Xj, QI)) +
                geom_line(aes(group = SimID, alpha = PercRank),
                    size = psize, colour = pcolour) +
                geom_smooth(method = method, size = lsize, se = FALSE,
                    color = lcolour) +
                geom_hline(aes(yintercept = 1), linetype = "dotted")
        }
    }
    # Plot ribbons
    else if (type == 'ribbons'){
        suppressWarnings(
        if (qi == "Hazard Rate"){
            if ('strata' %in% names(obj)) {
                obj <- MinMaxLines(df = obj, hr = TRUE, strata = TRUE)
                .e <- environment()
                p <- ggplot(obj, aes(x = Time, y = HRate,
                        colour = factor(HRValue), fill = factor(HRValue)),
                        environment = .e) +
                        geom_line(size = lsize) +
                        geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                            alpha = alpha, linetype = 0) +
                        geom_ribbon(aes(ymin = Min, ymax = Max),
                            alpha = alpha, linetype = 0) +
                        facet_grid(. ~ Strata) +
                        scale_colour_brewer(palette = spalette,
                            name = leg.name, guide = legend) +
                        scale_fill_brewer(palette = spalette,
                            name = leg.name, guide = legend)
            } else if (!('strata' %in% names(obj))){
                obj <- MinMaxLines(df = obj, hr = TRUE)
                .e <- environment()
                p <- ggplot(obj, aes(Time, Median, colour = factor(HRValue),
                              fill = factor(HRValue)), environment = .e) +
                    geom_line(size = lsize) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, linetype = 0) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                                    linetype = 0) +
                    scale_colour_brewer(palette = spalette,
                        name = leg.name) +
                    scale_fill_brewer(palette = spalette, name = leg.name)
            }
        } else if (qi == "First Difference"){
            obj <- MinMaxLines(df = obj)
            .e <- environment()
            p <- ggplot(obj, aes(Xj, Median), environment = .e) +
                geom_line(size = lsize, colour = lcolour) +
                geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                    alpha = alpha, fill = pcolour) +
                geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                                fill = pcolour) +
                geom_hline(aes(yintercept = 0), linetype = "dotted")
        } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
            obj <- MinMaxLines(df = obj)
            .e <- environment()
            p <- ggplot(obj, aes(Xj, Median), environment = .e) +
                geom_line(size = lsize, colour = lcolour) +
                geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                    alpha = alpha, fill = pcolour) +
                geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                                fill = pcolour) +
                geom_hline(aes(yintercept = 1), linetype = "dotted")
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
