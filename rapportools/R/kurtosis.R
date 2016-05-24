#' Kurtosis
#'
#' Calculates kurtosis coefficient for given variable (see \code{\link{is.variable}}), \code{matrix} or a \code{data.frame}.
#' @param x a \code{variable}, \code{matrix} or a \code{data.frame}
#' @param na.rm should \code{NA}s be removed before computation?
#' @references Tenjovic, L. (2000). Statistika u psihologiji - prirucnik. Centar za primenjenu psihologiju.
#' @examples
#' set.seed(0)
#' x <- rnorm(100)
#' kurtosis(x)
#' kurtosis(matrix(x, 10))
#' kurtosis(mtcars)
#' rm(x)
#' @export
kurtosis <- function(x, na.rm = TRUE){

    if (is.variable(x)){

        if (na.rm)
            x <- na.omit(x)

        m <- base::mean(x)
        s <- stats::sd(x)
        n <- length(x)
        (((base::sum((x - m) ^ 4) / n) / s ^ 4) - 3)
    } else {

        if (is.matrix(x))
            apply(x, 2, kurtosis, na.rm = na.rm)
        else if (is.data.frame(x))
            sapply(x, kurtosis, na.rm = na.rm)
        else
            stop('unsupported type')
    }
}


#' Skewness
#'
#' Calculates skewness coefficient for given variable (see \code{\link{is.variable}}), \code{matrix} or a \code{data.frame}.
#' @param x a \code{variable}, \code{matrix} or a \code{data.frame}
#' @param na.rm should \code{NA}s be removed before computation?
#' @references Tenjovic, L. (2000). Statistika u psihologiji - prirucnik. Centar za primenjenu psihologiju.
#' @examples
#' set.seed(0)
#' x <- rnorm(100)
#' skewness(x)
#' skewness(matrix(x, 10))
#' skewness(mtcars)
#' rm(x)
#' @export
skewness <- function(x, na.rm = TRUE){

    if (is.variable(x)){

        if (na.rm)
            x <- na.omit(x)

        m <- base::mean(x)
        s <- stats::sd(x)
        n <- length(x)
        (base::sum((x - m) ^ 3) / n) / s ^ 3

    } else {

        if (is.matrix(x))
            apply(x, 2, skewness, na.rm = na.rm)
        else if (is.data.frame(x))
            sapply(x, skewness, na.rm = na.rm)
        else
            stop('unsupported type')
    }
}
