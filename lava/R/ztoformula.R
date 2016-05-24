##' Converts strings to formula
##'
##' Converts a vector of predictors and a vector of responses (characters) i#nto
##' a formula expression.
##'
##'
##' @param y vector of predictors
##' @param x vector of responses
##' @return An object of class \code{formula}
##' @author Klaus K. Holst
##' @seealso \code{\link{as.formula}},
##' @keywords models utilities
##' @examples
##'
##' toformula(c("age","gender"), "weight")
##'
##' @export
toformula <- function (y = ".", x = ".")
{
    xst <- x[1]
    xn <- length(x)
    if (xn > 1)
        for (i in 2:length(x)) {
            xst <- paste(xst, "+", x[i])
        }
    yst <- y[1]
    yn <- length(y)
    if (yn > 1) {
        yst <- paste0("c(", yst)
        for (i in 2:length(y)) {
            yst <- paste0(yst, ", ", y[i])
        }
        yst <- paste0(yst, ")")
    }
    ff <- paste(yst, "~", xst)
    return(as.formula(ff))
}
