#' @title Simple Frequency Table
#'
#' @description Creates a simple frequency data.frame.
#'
#' @param x A vector of values for which the frequency is desired.
#' @param weighs A vector of weights.
#' @param breaks one of: 1) a vector giving the breakpoints between histogram
#'  cells; 2) a function to compute the vector of breakpoints; 3) a single
#'  number giving the number of cells for the histogram; 4) a character string
#'   naming an algorithm to compute the number of cells (see 'Details'); 5) a
#'    function to compute the number of cells.
#' @param digits The number of significant digits required.
#' @param include.lowest Logical; if \code{TRUE}, an x[i] equal to the breaks value will be included in the first (or last) category or bin.
#' @param order The order method.
#' @param useNA Logical; if \code{TRUE} NA's values are included.
#' @param \dots Additional arguements (currently ignored)
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#'
#' @seealso \code{\link{Frequency}}, \code{\link{crosstable}}.
#'
#' @examples
#' data(Presidents)
#'
#' freq(Presidents$winner.party)
#'
#'
#' @rdname freq
#' @export
`freq` <- function(x, weighs = NULL, breaks = graphics::hist(x, plot = FALSE)$breaks, digits=3, include.lowest = TRUE, order = c("desc", "asc","level", "name"), useNA = c("no", "ifany", "always"),...) UseMethod("freq")

#' @rdname freq
#' @export
`freq.default` <-
  function(x, weighs = NULL, breaks = graphics::hist(x, plot = FALSE)$breaks, digits=3, include.lowest = TRUE, order = c("desc", "asc","level", "name"), useNA = c("no", "ifany", "always"),...){

    # check if x is a vector (do not use is.vector())
    if(!(is.atomic(x) || is.list(x))) stop("'x' must be a vector")

    if(inherits(x, "table")){
      tab <- x

    } else {

      if(is.numeric(x)){
        x <- base::cut(x, breaks = breaks, include.lowest = include.lowest, ordered_result = TRUE,...)
      }

      tab <- base::table(x, useNA = useNA)
    }

    # how should the table be sorted, by name, level or frq? (NULL means "desc")
    switch(.Match(arg=order, choices = c( "desc", "asc", "level", "name")),
           level  = {  }
           , name   = { tab <- tab[rownames(tab)] }
           , asc    = { tab <- sort(tab) }
           , desc   = { tab <- -sort(-tab) }
    )

    ptab <- base::prop.table(tab)
    names(tab)[is.na(names(tab))] <- "<NA>"
    out <- data.frame(class = names(tab),
                      Freq = as.vector(tab[]), Prop = round(as.vector(ptab[]),digits))
    #cumfreq = cumsum(tab[]), cumperc = round(cumsum(ptab[]),digits))
    rownames(out) <- NULL # enumerate from 1:nrow(z)
    class(out) <- c("freq", "data.frame")
    return(out)
  }##-end of freq
NULL
