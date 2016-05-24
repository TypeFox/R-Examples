#' Plot simulated linear multiplicative interactions from Cox Proportional
#' Hazards Models
#'
#' \code{simGG.siminteract} uses ggplot2 to plot the quantities of
#' interest from \code{siminteract} objects, including marginal effects, first
#' differences, hazard ratios, and hazard rates.
#'
#' @param obj a \code{siminteract} class object
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
#' comparisons to plot. Not relevant for \code{qi = "Marginal Effect"}. Default
#' palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param legend specifies what type of legend to include (if applicable). The
#' default is \code{legend = "legend"}. To hide the legend use
#' \code{legend = FALSE}. See the \code{\link{discrete_scale}} for more details.
#' @param leg.name name of the legend (if applicable).
#' @param lcolour character string colour of the smoothing line. The default is
#' hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if
#' \code{qi = "Marginal Effect"}.
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
#'  \code{\link{coxsimSpline}}) as well as the central 50 percent of this area.
#' It also plots a line for the median value of the full area, so values in
#' \code{method} are ignored. One of the key advantages of using ribbons
#' rather than points is that it creates plots with smaller file sizes.
#' @param ... Additional arguments. (Currently ignored.)
#'
#' @return a \code{gg} \code{ggplot} class object
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
#' # Simulate Marginal Effect of lethal for multiple values of prevgenx
#' Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx",
#'                      X2 = seq(2, 115, by = 2), nsim = 100)
#'
#' # Plot quantities of interest
#' simGG(Sim1)
#' simGG(Sim1, rug_position = 'jitter')
#'
#' \dontrun{
#' # Change the order of the covariates to make a more easily
#' # interpretable hazard ratio graph.
#' M2 <- coxph(Surv(acttime, censor) ~ prevgenx*lethal,
#'              data = CarpenterFdaData)
#'
#' # Simulate Hazard Ratio of lethal for multiple values of prevgenx
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
#'                        X1 = 100, X2 = 1, qi = "Hazard Rate")
#'
#' # Plot quantities of interest
#' simGG(Sim1, xlab = "\nprevgenx",
#'        ylab = "Marginal Effect of lethal\n")
#' simGG(Sim2, type = 'ribbons', rug_position = 'jitter')
#' simGG(Sim3)
#' simGG(Sim4, to = 150, type = 'lines', legend = FALSE)
#' }
#'
#' @details Uses ggplot2 to plot the quantities of interest from
#' \code{siminteract} objects, including marginal effects, first differences,
#' hazard ratios, and hazard rates. If there are multiple strata, the quantities
#' of interest will be plotted in a grid by strata.
#'
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no
#' effect, for time-varying hazard ratio graphs. No line is created for hazard
#' rates.
#'
#' Note: if \code{qi = "Hazard Ratio"} or \code{qi = "First Difference"} then
#' you need to have choosen more than one fitted value for \code{X1} in
#' \code{\link{coxsimInteract}}.
#'
#'
#' @seealso \code{\link{coxsimInteract}}, \code{\link{simGG.simlinear}},
#' and \code{ggplot2}
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006.
#' ''Understanding Interaction Models: Improving Empirical Analyses.'' Political
#'  Analysis 14(1): 63-82.
#'
#' Keele, Luke. 2010. ''Proportionally Difficult: Testing for Nonproportional
#' Hazards in Cox Models.'' Political Analysis 18(2): 189-205.
#'
#' Carpenter, Daniel P. 2002. ''Groups, the Media, Agency Waiting Costs, and
#' FDA Drug Approval.'' American Journal of Political Science 46(3): 490-505.
#'
#' @import ggplot2
#' @import mgcv
#' @importFrom dplyr rename
#'
#' @method simGG siminteract
#' @export

simGG.siminteract <- function(obj, from = NULL,
                            rug = TRUE, rug_position = "identity",
                            to = NULL, xlab = NULL,
                            ylab = NULL, title = NULL, method = "auto",
                            spalette = "Set1", legend = "legend",
                            leg.name = "", lcolour = "#2B8CBE", lsize = 1,
                            pcolour = "#A6CEE3", psize = 1, alpha = 0.2,
                            type = "ribbons", ...)
    {
    HRate <- Lower50 <- Upper50 <- Min <- Max <- Median <-
    Time <- QI <- HRValue <- X1 <- X2 <- SimID <- xaxis <- NULL
    if (!inherits(obj, "siminteract")){
        stop("must be a siminteract object", call. = FALSE)
    }
    if (type == 'ribbons' & method != "auto"){
        message("The method argument is ignored if ribbons = TRUE. Central tendency summarised with the median.")
    }
    # Find quantity of interest
    qi <- class(obj)[[2]]

    # Create y-axis label
    if (is.null(ylab)) ylab <- paste(qi, "\n")

    # Create x-axis label
    if (qi != "Hazard Rate"){
        if (qi == "Marginal Effect") {
            b_for_x <- 'b2'
            x_for_x <- 'X2'
        }
        else if (qi == "First Difference" | qi == "Hazard Ratio"){
            b_for_x <- 'b1'
            x_for_x <- 'X1'
        }
        x_var <- attr(obj, b_for_x)
        if (is.null(xlab)) xlab <- paste("\n", x_var)

        # Extract rug values
        rugger <- rugExtract(obj, x = x_for_x, rug_var = x_var)
    }

    # Convert obj to data frame
    obj <- as.data.frame(obj)

    # Drop simulations that include outliers
    if (type == 'lines'){
        obj <- OutlierDrop(obj)
    }

    # Alpha gradient based on percentile in the distribution
    if (type != 'ribbons' & qi != 'Hazard Rate' & qi != 'Marginal Effect'){
        obj <- PercRank(obj, xaxis = 'X1')
    } else if (type != 'ribbons' & qi == 'Marginal Effect'){
        obj <- PercRank(obj, xaxis = 'X2')
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
                obj$HRValue <- as.factor(obj$HRValue)
                .e <- environment()
                p <- ggplot(obj, aes(x = Time, y = HRate,
                            colour = HRValue), environment = .e) +
                    geom_point(aes(alpha = PercRank), size = psize) +
                    geom_smooth(method = method, size = lsize, se = FALSE) +
                    facet_grid(.~ Strata) +
                    scale_colour_brewer(palette = spalette, name = leg.name,
                                        guide = legend) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
            } else if (!('strata' %in% names(obj))){
                obj$HRValue <- as.factor(obj$HRValue)
                .e <- environment()
                p <- ggplot(obj, aes(x = Time, y = HRate,
                              colour = HRValue), environment = .e) +
                        geom_point(shape = 21, aes(alpha = PercRank),
                            size = psize) +
                        geom_smooth(method = method, size = lsize,
                            se = FALSE) +
                        scale_colour_brewer(palette = spalette,
                            name = leg.name, guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        }
        else if (qi == "Marginal Effect"){
            p <- ggplot(obj, aes(X2, QI)) +
                    geom_point(shape = 21, aes(alpha = PercRank), size = psize,
                               colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                                color = lcolour) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        }
        else if (qi == "First Difference"){
            X1Unique <- obj[!duplicated(obj[, "X1"]), ]
            if (nrow(X1Unique) <= 1){
                message("X1 must have more than one fitted value.")
            } else {
                obj$X2 <- as.factor(obj$X2)
                .e <- environment()
                p <- ggplot(obj, aes(X1, QI, colour = X2, group = X2),
                            environment = .e) +
                        geom_point(shape = 21, aes(alpha = PercRank),
                                                  size = psize) +
                        geom_smooth(method = method, size = lsize, se = FALSE) +
                        geom_hline(aes(yintercept = 0), linetype = "dotted") +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                                            guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        }
        else if (qi == "Hazard Ratio"){
            X1Unique <- obj[!duplicated(obj[, "X1"]), ]
            if (nrow(X1Unique) <= 1){
                message("X1 must have more than one fitted value.")
            } else {
                obj$X2 <- as.factor(obj$X2)
                .e <- environment()
                p <- ggplot(obj, aes(X1, QI, colour = X2, group = X2),
                        environment = .e) +
                        geom_point(shape = 21, aes(alpha = PercRank),
                           size = psize) +
                        geom_smooth(method = method, size = lsize, se = FALSE) +
                        geom_hline(aes(yintercept = 1), linetype = "dotted") +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                            guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        }
    }
    # Plot lines
    else if (type == 'lines'){
        if (qi == "Hazard Rate"){
            if ('strata' %in% names(obj)) {
                obj$HRValue <- as.factor(obj$HRValue)
                .e <- environment()
                p <- ggplot(obj, aes(x = Time, y = HRate, colour = HRValue),
                        environment = .e) +
                        geom_line(aes(group = interaction(SimID, HRValue),
                                      alpha = PercRank), size = psize) +
                        geom_smooth(aes(group = HRValue),
                            method = method, size = lsize, se = FALSE) +
                        facet_grid(.~ Strata) +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                                            guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            } else if (!('strata' %in% names(obj))){
                obj$HRValue <- as.factor(obj$HRValue)
                .e <- environment()
                p <- ggplot(obj, aes(x = Time, y = HRate, colour = HRValue),
                        environment = .e) +
                        geom_line(aes(group = SimID, alpha = PercRank),
                            shape = 21, size = psize) +
                        geom_smooth(aes(group = HRValue),
                            method = method, size = lsize, se = FALSE) +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                            guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        }
        else if (qi == "Marginal Effect"){
            p <- ggplot(obj, aes(X2, QI)) +
                    geom_line(aes(group = SimID, alpha = PercRank),
                        size = psize, colour = pcolour) +
                    geom_smooth(method = method, size = lsize, se = FALSE,
                                color = lcolour) +
                    scale_alpha_continuous(range = c(0, alpha), guide = FALSE)
        }
        else if (qi == "First Difference"){
            X1Unique <- obj[!duplicated(obj[, "X1"]), ]
            if (nrow(X1Unique) <= 1){
                message("X1 must have more than one fitted value.")
            } else {
                obj$X2 <- as.factor(obj$X2)
                .e <- environment()
                p <- ggplot(obj, aes(X1, QI, colour = X2), environment = .e) +
                        geom_line(aes(group = interaction(SimID, X2),
                                      alpha = PercRank), size = psize) +
                        geom_smooth(aes(group = X2), method = method,
                                        size = lsize, se = FALSE) +
                        geom_hline(aes(yintercept = 0), linetype = "dotted") +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                                            guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        }
        else if (qi == "Hazard Ratio"){
            X1Unique <- obj[!duplicated(obj[, "X1"]), ]
            if (nrow(X1Unique) <= 1){
                message("X1 must have more than one fitted value.")
            } else {
                obj$X2 <- as.factor(obj$X2)
                p <- ggplot(obj, aes(X1, QI, colour = X2)) +
                        geom_line(aes(group = interaction(SimID, X2),
                                      alpha = PercRank), size = psize) +
                        geom_smooth(aes(group = X2), method = method,
                                        size = lsize, se = FALSE) +
                        geom_hline(aes(yintercept = 1), linetype = "dotted") +
                        scale_colour_brewer(palette = spalette, name = leg.name,
                                            guide = legend) +
                        scale_alpha_continuous(range = c(0, alpha),
                            guide = FALSE)
            }
        }
    }
    # Plot ribbons
    else if (type == 'ribbons'){
        suppressWarnings(
        if (qi == "Hazard Rate"){
            if ('strata' %in% names(obj)) {
                obj <- MinMaxLines(df = obj, hr = TRUE, strata = TRUE)
                obj$HRValue <- as.factor(obj$HRValue)
                .e <- environment()
                p <- ggplot(obj, aes(x = Time, y = HRate,
                    colour = HRValue, fill = HRValue),
                    environment = .e) +
                        geom_line(size = lsize, alpha = alpha) +
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
                obj$HRValue <- as.factor(obj$HRValue)
                .e <- environment()
                p <- ggplot(obj, aes(Time, Median, colour = HRValue,
                        fill = HRValue), environment = .e) +
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
        }
        else if (qi == "Marginal Effect"){
            obj <- MinMaxLines(df = obj, byVars = c("X2"))
            .e <- environment()
            p <- ggplot(obj, aes(X2, Median), environment = .e) +
                    geom_line(size = lsize, colour = lcolour) +
                    geom_ribbon(aes(ymin = Lower50, ymax = Upper50),
                        alpha = alpha, fill = pcolour) +
                    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha,
                                    fill = pcolour)
        }
        else if (qi == "Hazard Ratio" | qi == "First Difference"){
            X1Unique <- obj[!duplicated(obj[, "X1"]), ]
            if (nrow(X1Unique) <= 1){
                message("X1 must have more than one fitted value.")
            } else {
                obj <- MinMaxLines(df = obj, byVars = c("X1", "X2"))
                obj$X2 <- as.factor(obj$X2)
                .e <- environment()
                p <- ggplot(obj, aes(X1, Median, colour = X2,
                            fill = X2), environment = .e) +
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
            }
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
