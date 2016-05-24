#' @title plot the result of univariate regression result
#'
#' @description plot coefficients, OR or HR of univariate regressions.

#' @param x A reg object
#' @param limits A numeric vector of length two providing limits of the scale. Use NA to refer to the existing minimum or maximum value.
#' @param sort A character determining the order of variables  to plot,'alphabetical' or 'order'. The later is the default to sort varibales according to their effect size.
#' @param \dots additional arguments
#' @seealso \code{\link{reg}}, \code{\link{dataframe.reg}},  \code{\link{detail.reg}}
#' @import ggplot2
#' @export
#' @examples
#' reg_glm<-reg(data = diabetes, y = 5, factor = c(1, 3, 4), model = 'glm')
#' plot(reg_glm)
#' plot(reg_glm, limits = c(NA, 3))


plot.reg <- function(x, limits = c(NA, NA), sort = "order", ...) {
    if (class(x) != "reg") {
        stop("x should be a `reg` object.", call. = FALSE)
    }
    stopifnot((is.na(limits[1]) || is.numeric(limits[1])) && (is.na(limits[2]) ||
        is.numeric(limits[2])))

    data <- x$dataframe
    names(data)[(NCOL(x$dataframe) - 2):NCOL(x$dataframe)] <- c("center", "low",
        "high")

    if (is.na(limits[1]))
        limits[1] <- min(data$low, na.rm = TRUE)
    if (is.na(limits[2]))
        limits[2] <- max(data$high, na.rm = TRUE)

    if (!is.na(diff(limits)) && diff(limits) < 0) {
        stop("A problem with limits, please check it out!", call. = FALSE)
    }
    adj <- diff(limits) * 0.05

    data <- cbind(data, beyond1 = ifelse(data$low < limits[1] - adj, TRUE,
        NA), beyond2 = ifelse(data$high > limits[2] + adj, TRUE, NA))

    data$low <- ifelse(data$low < limits[1], limits[1], data$low)
    data$high <- ifelse(data$high > limits[2], limits[2], data$high)

    if (sort == "order") {
        data <- data[order(data$center, decreasing = TRUE), ]
        data$term <- factor(data$term, levels = data$term)
    } else if (sort == "alphabetical") {
        data <- data[order(data$term, decreasing = TRUE), ]
        data$term <- factor(data$term, levels = data$term)
    }

    p <- ggplot(data, aes(term, center, ymin = low, ymax = high, color = term)) +
        geom_linerange() + geom_point(size = 1.5, shape = 15) + theme(plot.title = element_text(size = 30),
        axis.text.x = element_text(hjust = 1, vjust = 0, size = 8), legend.position = "none") +
        coord_flip(ylim = limits) + geom_segment(aes(x = term, y = low * beyond2,
        xend = term, yend = limits[2] * beyond2 + adj, group = term), arrow = arrow(type = "closed",
        length = unit(0.15, "inches")), na.rm = TRUE) + geom_segment(aes(x = term,
        y = high * beyond1, xend = term, yend = limits[1] * beyond1 - adj,
        group = term), arrow = arrow(type = "closed", length = unit(0.15, "inches")),
        na.rm = TRUE)

    p <- switch(x$detail$call$model,
                lm = p + ylab("coefficients(95%CI)")+geom_abline(intercept = 0, slope = 0, colour = "red", size = 0.1),
                glm = p + ylab("OR(95%CI)")+geom_abline(intercept = 1, slope = 0, colour = "red", size = 0.1),
                coxph = p + ylab("HR(95%CI)")+geom_abline(intercept = 1, slope = 0, colour = "red", size = 0.1))
    return(p)
}


utils::globalVariables(c("term", "center", "low", "high", "beyond1", "beyond2"))
