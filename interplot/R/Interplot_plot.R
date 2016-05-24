#' Plot Conditional Coefficients in Models with Interaction Terms
#' 
#' Graph based on the data frame of statistics about the conditional effect of an interaciton.
#' 
#' @param m A model object including an interaction term, or, alternately, a data frame recording conditional coefficients. This data frame should includes four columns:
#' \itemize{
#'    \item fake: The sequence of \code{var1} (the item whose effect will be conditioned on in the interaction);
#'    \item coef1: The point estimates of the coefficient of \code{var1} at each break point.
#'    \item ub: The upper bound of the simulated 95\% CI.
#'    \item lb: The lower bound of the simulated 95\% CI.
#' }
#' @param var1 The name (as a string) of the variable of interest in the interaction term; its conditional coefficient estimates will be plotted.
#' @param var2 The name (as a string) of the other variable in the interaction term.
#' @param plot A logical value indicating whether the output is a plot or a dataframe including the conditional coefficient estimates of var1, their upper and lower bounds, and the corresponding values of var2.
#' @param hist A logical value indicating if there is a histogram of `var2` added at the bottom of the conditional effect plot.
#' @param var2_dt A numerical value indicating the frequency distibution of `var2`. It is only used when `hist == TRUE`. When the object is a model, the default is the distribution of `var2` of the model. 
#' @param point A logical value determining the format of plot. By default, the function produces a line plot when var2 takes on ten or more distinct values and a point (dot-and-whisker) plot otherwise; option TRUE forces a point plot.
#' @param sims Number of independent simulation draws used to calculate upper and lower bounds of coefficient estimates: lower values run faster; higher values produce smoother curves.
#' @param xmin A numerical value indicating the minimum value shown of x shown in the graph. Rarely used.
#' @param xmax A numerical value indicating the maximum value shown of x shown in the graph. Rarely used.
#' @param ercolor A character value indicating the outline color of the whisker or ribbon.
#' @param esize A numerical value indicating the size of the whisker or ribbon.
#' @param ralpha A numerical value indicating the transparency of the ribbon.
#' @param rfill A character value indicating the filling color of the ribbon.
#' @param ... Other ggplot aesthetics arguments for points in the dot-whisker plot or lines in the line-ribbon plots. Not currently used.
#' 
#' @details \code{interplot.plot} is a S3 method from the \code{interplot}. It generates plots of conditional coefficients.
#' 
#' Because the output function is based on \code{\link[ggplot2]{ggplot}}, any additional arguments and layers supported by \code{ggplot2} can be added with the \code{+}. 
#' 
#' @return The function returns a \code{ggplot} object.
#' 
#' @import  ggplot2
#' @importFrom graphics hist
#' 
#' @export

## S3 method for class 'data.frame'
interplot.plot <- function(m, var1, var2, plot = TRUE, hist = FALSE, var2_dt = NULL, point = FALSE, 
    sims = 5000, xmin = NA, xmax = NA, ercolor = NA, esize = 0.5, ralpha = 0.5, rfill = "grey70", 
    ...) {
    steps <- nrow(m)
    levels <- sort(unique(m$fake))
    ymin <- ymax <- vector() # to deal with the "no visible binding for global variable" issue
    if (hist == FALSE) {
        if (steps < 10 | point == T) {
            if (is.na(ercolor)) 
                {
                  ercolor <- "black"
                }  # ensure whisker can be drawn
            coef.plot <- ggplot(m, aes_string(x = "fake", y = "coef1")) + geom_point(...) + geom_errorbar(aes_string(ymin = "lb", 
                ymax = "ub"), width = 0, color = ercolor, size = esize) + scale_x_continuous(breaks = levels) + 
                ylab(NULL) + xlab(NULL)
        } else {
            coef.plot <- ggplot(m, aes_string(x = "fake", y = "coef1")) + geom_line(...) + geom_ribbon(aes_string(ymin = "lb", 
                ymax = "ub"), alpha = ralpha, color = ercolor, fill = rfill) + ylab(NULL) + xlab(NULL)
        }
        return(coef.plot)
    } else {
        if (steps < 10 | point == T) {
            if (is.na(ercolor)) 
                {
                  ercolor <- "black"
                }  # ensure whisker can be drawn
            
            yrange <- c(m$ub, m$lb, var2_dt)
            maxdiff <- (max(yrange) - min(yrange))
            
            break_var2 <- length(unique(var2_dt))
            if (break_var2 >= 100) 
                break_var2 <- 100
            hist.out <- hist(var2_dt, breaks = break_var2, plot = FALSE)
            
            n.hist <- length(hist.out$mids)
            dist <- hist.out$mids[2] - hist.out$mids[1]
            hist.max <- max(hist.out$counts)
            
            histX <- data.frame(ymin = rep(min(yrange) - maxdiff/5, n.hist), ymax = hist.out$counts/hist.max * 
                maxdiff/5 + min(yrange) - maxdiff/5, xmin = hist.out$mids - dist/2, xmax = hist.out$mids + 
                dist/2)
            
            coef.plot <- ggplot()
            coef.plot <- coef.plot + geom_rect(data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, 
                ymax = ymax), colour = "gray50", alpha = 0, size = 0.5)  #histgram
            
            coef.plot <- coef.plot + geom_point(data = m, aes_string(x = "fake", y = "coef1")) + 
                geom_errorbar(data = m, aes_string(x = "fake", ymin = "lb", ymax = "ub"), width = 0, 
                  color = ercolor, size = esize) + scale_x_continuous(breaks = levels) + ylab(NULL) + 
                xlab(NULL)
            
        } else {
            
            yrange <- c(m$ub, m$lb)
            
            maxdiff <- (max(yrange) - min(yrange))
            
            break_var2 <- length(unique(var2_dt))
            if (break_var2 >= 100) 
                break_var2 <- 100
            hist.out <- hist(var2_dt, breaks = break_var2, plot = FALSE)
            
            n.hist <- length(hist.out$mids)
            dist <- hist.out$mids[2] - hist.out$mids[1]
            hist.max <- max(hist.out$counts)
            
            histX <- data.frame(ymin = rep(min(yrange) - maxdiff/5, n.hist), ymax = hist.out$counts/hist.max * 
                maxdiff/5 + min(yrange) - maxdiff/5, xmin = hist.out$mids - dist/2, xmax = hist.out$mids + 
                dist/2)
            
            
            
            # interplot.plot(m = coef, var1 = 'cyl', var2 = 'wt')
            
            coef.plot <- ggplot()
            coef.plot <- coef.plot + geom_rect(data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, 
                ymax = ymax), colour = "gray50", alpha = 0, size = 0.5)
            
            
            coef.plot <- coef.plot + geom_line(data = m, aes_string(x = "fake", y = "coef1")) + 
                geom_ribbon(data = m, aes_string(x = "fake", ymin = "lb", ymax = "ub"), alpha = ralpha, 
                  color = ercolor, fill = rfill) + ylab(NULL) + xlab(NULL)
        }
        return(coef.plot)
    }
    
} 
