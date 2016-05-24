##
#' Generic method to make a footnote indicating significant values.
#'
#' @param data Dataset with statistics.
#' @return \code{apa.signif} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{signif}{\code{pot {ReporteRs}} object}
#' @importFrom "ReporteRs" "textProperties" "pot"
#' @export
#'
#' @examples
#'
#' # Specify statistics
#' example <- data.frame(
#'   c("Column 1", "Column 2", "Column 3"),
#'   c(3.45, 5.21, 2.64),
#'   c("**", "", "***")
#' )
#'
#' # Use apa.descriptives function
#' apa.signif(data = example)
##
apa.signif = function(data=data.frame()) UseMethod("apa.signif")

##
#' Default method to make a footnote indicating significant values.
#'
#' @param data Dataset with statistics.
#' @return \code{apa.signif} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{signif}{\code{pot {ReporteRs}} object}
#' @importFrom "ReporteRs" "textProperties" "pot"
#' @export
#'
#' @examples
#'
#' # Specify statistics
#' example <- data.frame(
#'   c("Column 1", "Column 2", "Column 3"),
#'   c(3.45, 5.21, 2.64),
#'   c("**", "", "***")
#' )
#'
#' # Use apa.descriptives function
#' apa.signif(data = example)
##
apa.signif.default = function(data=data.frame()) {

  est = apaSignificance(data)
  est$call = match.call()
  class(est) = "apa.signif"
  est

}

##
#' Define a print method
#'
#' @param  x A \code{apa.signif} object
#' @export
##
print.apa.signif = function(x, ...) {
  if(x$succes == TRUE) {
    cat("\n")
    cat("Succesfully generated significance footnote.")
    cat("\n\n")
  }
}

# The main function
apaSignificance = function(data) {

  # Initialize function
  options(warn = 0)

  # Check if a valid data frame is supplied
  if ((!is.data.frame(data)) || (is.data.frame(data) && nrow(data) == 0)) {
    error = "Invalid data is supplied."
    warning(error)
    return(list(succes = error))
  }

  # Check the size of the dataset
  if (ncol(data) > 20 | nrow(data) > 100) {
    error = "The supplied data is too big to generate an APA formatted table."
    warning(error)
    return(list(succes = error))
  } else {

    # Convert factors to characters
    i = sapply(data, is.factor)
    data[i] = lapply(data[i], as.character)

    # Convert "+" symbol to unicode dagger symbol
    data[which(data == "+", arr.ind = TRUE)] = "\u2020"

    has.signif1 = apply(data, c(1, 2), function(x) any(x == "\u2020"))
    has.signif2 = apply(data, c(1, 2), function(x) any(x == "*"))
    has.signif3 = apply(data, c(1, 2), function(x) any(x == "**"))
    has.signif4 = apply(data, c(1, 2), function(x) any(x == "***"))

    if (TRUE %in% has.signif1) {
      sig1 = ReporteRs::pot("\u2020", ReporteRs::textProperties(font.family = "Times", font.size = 12, vertical.align = "superscript")) +
        ReporteRs::pot("p", ReporteRs::textProperties(font.family = "Times", font.size = 12, font.style = "italic")) +
        ReporteRs::pot(" < .10", ReporteRs::textProperties(font.family = "Times", font.size = 12))
    } else {
      sig1 = ReporteRs::pot("")
    }

    if ((TRUE %in% has.signif2) || (TRUE %in% has.signif3 || TRUE %in% has.signif4)) {
      if(!"" %in% sig1[[1]]$value) { sig1 = sig1 + ReporteRs::pot("; ", ReporteRs::textProperties(font.family = "Times", font.size = 12)) }
      sig2 = ReporteRs::pot("*", ReporteRs::textProperties(font.family = "Times", font.size = 12)) +
        ReporteRs::pot("p", ReporteRs::textProperties(font.family = "Times", font.size = 12, font.style = "italic")) +
        ReporteRs::pot(" < .05", ReporteRs::textProperties(font.family = "Times", font.size = 12))
    } else {
      sig2 = ReporteRs::pot("")
    }

    if ((TRUE %in% has.signif3) || (TRUE %in% has.signif4)) {
      if(!"" %in% sig2[[1]]$value) { sig2 = sig2 + ReporteRs::pot("; ", ReporteRs::textProperties(font.family = "Times", font.size = 12)) }
      sig3 = ReporteRs::pot("**", ReporteRs::textProperties(font.family = "Times", font.size = 12)) +
        ReporteRs::pot("p", ReporteRs::textProperties(font.family = "Times", font.size = 12, font.style = "italic")) +
        ReporteRs::pot(" < .01", ReporteRs::textProperties(font.family = "Times", font.size = 12))
    } else {
      sig3 = ReporteRs::pot("")
    }

    if (TRUE %in% has.signif4) {
      if(!"" %in% sig3[[1]]$value) { sig3 = sig3 + ReporteRs::pot("; ", ReporteRs::textProperties(font.family = "Times", font.size = 12)) }
      sig4 = ReporteRs::pot("***", ReporteRs::textProperties(font.family = "Times", font.size = 12)) +
        ReporteRs::pot("p", ReporteRs::textProperties(font.family = "Times", font.size = 12, font.style = "italic")) +
        ReporteRs::pot(" < .001", ReporteRs::textProperties(font.family = "Times", font.size = 12))
    } else {
      sig4 = ReporteRs::pot("")
    }
    apa.signif =  sig1 + sig2 + sig3 + sig4

  }

  return(list(succes = TRUE, signif = apa.signif))

}
