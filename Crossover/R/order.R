#' Sorts sequences of a design into a canonical order
#' 
#' Sorts sequences of a design into a canonical order.
#' 
#' When comparing bigger designs this ordering easily allows to check whether two
#' designs are equal.
#' 
#' @param design Cross-over design.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @keywords misc
#' @examples
#' 
#' getDesign("switchback5t")
#' canonicalOrder(getDesign("switchback5t"))
#' 
#' @export canonicalOrder
canonicalOrder <- function(design) {
  design[, do.call(order, lapply(1:NROW(design), function(i) design[i, ]))]
}