##' format p values
##'
##' @param p p values
##' @param digits number of digits
##' @return formated p values
##' @export
##' @author David Hajage
plim <- function(p, digits = 4) {
  pround <- round(p, digits)
  lim <- 10^(-digits)
  ptxt <- vector("character", length(p))
  ptxt[pround < lim] <- paste("<", "0.", paste(rep("0", digits-1), collapse = ""), "1", sep = "")
  ptxt[pround >= lim] <- formatC(pround[pround >= lim], format = "f", digits = digits)
  return(ptxt)
}
