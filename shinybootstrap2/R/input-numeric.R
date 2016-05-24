#' Create a numeric input control
#'
#' Create an input control for entry of numeric values
#'
#' @param inputId Input variable to assign the control's value to
#' @param label Display label for the control
#' @param value Initial value
#' @param min Minimum allowed value
#' @param max Maximum allowed value
#' @param step Interval to use when stepping between min and max
#' @return A numeric input control that can be added to a UI definition.
#'
#' @family input elements
#' @seealso \code{\link[shiny]{updateNumericInput}}
#'
#' @examples
#' numericInput("obs", "Observations:", 10,
#'              min = 1, max = 100)
#' @export
numericInput <- function(inputId, label, value, min = NA, max = NA, step = NA) {

  # build input tag
  inputTag <- tags$input(id = inputId, type = "number", value = formatNoSci(value))
  if (!is.na(min))
    inputTag$attribs$min = min
  if (!is.na(max))
    inputTag$attribs$max = max
  if (!is.na(step))
    inputTag$attribs$step = step

  tagList(
    label %AND% tags$label(label, `for` = inputId),
    inputTag
  )
}
