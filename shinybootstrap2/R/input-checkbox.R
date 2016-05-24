#' Checkbox Input Control
#'
#' Create a checkbox that can be used to specify logical values.
#'
#' @param inputId Input variable to assign the control's value to.
#' @param label Display label for the control.
#' @param value Initial value (\code{TRUE} or \code{FALSE}).
#' @return A checkbox control that can be added to a UI definition.
#'
#' @family input elements
#' @seealso \code{\link{checkboxGroupInput}}, \code{\link[shiny]{updateCheckboxInput}}
#'
#' @examples
#' checkboxInput("outliers", "Show outliers", FALSE)
#' @export
checkboxInput <- function(inputId, label, value = FALSE) {
  inputTag <- tags$input(id = inputId, type="checkbox")
  if (!is.null(value) && value)
    inputTag$attribs$checked <- "checked"
  tags$label(class = "checkbox", `for` = inputId, inputTag, tags$span(label))
}
