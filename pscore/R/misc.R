#' Density Plot for a Long Dataset
#'
#' Internal function only, not meant for general use
#' Simple wrapper around ggplot2 functionaly to create
#' density plots, potentially for many variables and coloured by group.
#'
#' @param data A dataset (or melt()ed dataset)
#' @param melt Logical whether to melt() dataset
#' @param x name of variable for density
#' @param facet A variable to use for facetting
#' @param g A variable to use for grouping/colouring.  If \code{melt=TRUE}, this is
#'   used as id.var as well.
#' @param hist Logical whether to make a density plot or histogram (if TRUE).
#' @return A ggplot2 graph.
#' @import ggplot2 reshape2
#' @examples
#' # simple facetted plot
#' pscore:::ldensity(mtcars, TRUE)
#' # simple coloured plot
#' pscore:::ldensity(mtcars, x = "mpg", g = "factor(cyl)")
ldensity <- function(data, melt = FALSE, x, facet, g, hist=FALSE) {
    data <- as.data.frame(data)
    if (melt) {
        if (!missing(g)) {
            data <- melt(data, id.var = g)
        } else {
            data <- melt(data)
        }
        x <- "value"
        facet <- ~variable
    }

    if (!missing(g)) {
        p <- ggplot(data, aes_string(x = x, color = g))
    } else {
        p <- ggplot(data, aes_string(x = x))
    }

    if (hist) {
        p <- p + geom_histogram()
    } else {
        p <- p + geom_density()
    }

    if (!missing(facet)) {
        p <- p + facet_wrap(facet, scales = "free")
    }
    p + theme_bw()
}

#' Winsorize at specified percentiles
#'
#' Simple function winsorizes data at the specified percentile.
#'
#' @param d A vector, matrix, or data frame to be winsorized
#' @param percentile The percentile bounded by [0, 1] to winsorize data at
#' @param na.rm A logical whether to remove NAs.
#' @return winsorized data. Attributes are included to list the exact values
#'   (for each variable, if a data frame or matrix) used to winsorize
#'   at the lower and upper ends.
#' @export
#' @examples
#' dev.new(width = 10, height = 5)
#' par(mfrow = c(1, 2))
#' hist(as.vector(eurodist), main = "Eurodist")
#' hist(winsorizor(as.vector(eurodist), .05), main = "Eurodist with lower and upper\n5% winsorized")
winsorizor <- function(d, percentile, na.rm = TRUE) {
    stopifnot(percentile >= 0 && percentile <= 1)
    stopifnot(is.vector(d) || is.matrix(d) || is.data.frame(d))

    f <- function(x, percentile, na.rm) {
          low <- quantile(x, probs = 0 + percentile, na.rm = na.rm)
          high <- quantile(x, probs = 1 - percentile, na.rm = na.rm)

          out <- pmin(pmax(x, low), high)

          new.attr <- data.frame(low = low, high = high, percentile = percentile)
          rownames(new.attr) <- NULL

          attributes(out) <- c(attributes(x), winsorized = list(new.attr))

          return(out)
    }

    if (is.vector(d)) {
        out <- f(d, percentile = percentile, na.rm = na.rm)
    } else if (is.matrix(d) || is.data.frame(d)) {

        tmp <- lapply(1:ncol(d), function(i) f(d[, i], percentile = percentile, na.rm = na.rm))
        all.attr <- do.call(rbind, lapply(tmp, function(x) attr(x, "winsorized")))
        all.attr$variable <- colnames(d)
        rownames(all.attr) <- NULL

        if (is.matrix(d)) {
            out <- as.matrix(as.data.frame(tmp))
        } else {
            out <- as.data.frame(tmp)
        }
        attributes(out) <- c(attributes(d), winsorized = list(all.attr))
    }

    return(out)
}


#' @name BioDB
#' @title Biomarker Thresholds Database
#' @description This data set lists the clinical high risk thresholds for a variety of biomarkers
#' @docType data
#' @format a \code{list}.
#' @source Various publications
NULL
## BioDB <- list(
##     MetabolicSyndrome = data.frame(
##         Biomarker = c("SBP", "DBP", "Triglycerides", "HDL-C", "WaistCircumference", "BMI", "Glucose", "HbA1c"),
##         Units = c("mmHg", "mmHg", "mmol/L", "mmol/L", "cm", "kg/m^2", "mmol/L", "percent"),
##         Female = c(130, 85, 1.7, 1.3, 80, 25, 5.6, 5.7),
##         Male = c(130, 85, 1.7, 1.0, 94, 25, 5.6, 5.7),
##         Higherisbetter = c(0, 0, 0, 1, 0, 0, 0, 0)))
## save(BioDB, file = "../data/BioDB.rda")
