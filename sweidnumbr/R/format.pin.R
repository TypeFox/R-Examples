#' Formatting pin
#' 
#' Format \code{pin} for pretty printing
#'
#' @param x vector of class "pin" (see \code{\link{as.pin}}) or a vector that can be coerced to such
#' @param format. character string specifying the output format. 
#' \code{\%N} is used as a reference for the last four digits of the pin.
#' Format of the date is handled via \code{\link{strptime}}.
#' (\code{"\%Y\%m\%d\%N"} by default). \code{\%P} is an available  
#' shorthand for \code{"(\%C) \%y-\%m-\%d - \%N"}, a format aimed for 
#' maximal readability when used in long lists
#' @param ... arguments passed to \code{\link{format.Date}}
#'
#' @return character vector of same length as \code{x}
#' @export
#'
#' @examples
#' x <- as.pin(fake_pins$pin[1:10])
#' 
#' # Separate elements with hyphens:
#' format_pin(x, "%Y-%m-%d-%N")
#' 
#' 
#' # Separate even further
#' format_pin(x, "%C-%y-%m-%d-%N")
#'
#' # The special P-format for maximal readability
#' format_pin(x, "%P") 
#' 
#' # A custom representation
#' format_pin(x, "Borned %d of %B in %Y (a %A in week %U) with suffix no: %N")
#' 
#' # Extract only the year
#' format_pin(x, "%Y")
format_pin <- function(x, format. = "%Y%m%d%N", ...) {
  if (format. == "%P") format. <- "(%C) %y-%m-%d - %N"
  gsub_v <- Vectorize(gsub, "replacement")
  f <- gsub_v("%N", substr(x, 9, 12), format.)
  mapply(format, pin_to_date(x), f)
}
