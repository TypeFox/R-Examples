#' Create a text input control
#'
#' Create an input control for entry of unstructured text values
#'
#' @param inputId Input variable to assign the control's value to
#' @param label Display label for the control
#' @param value Initial value
#' @return A text input control that can be added to a UI definition.
#'
#' @family input elements
#' @seealso \code{\link[shiny]{updateTextInput}}
#'
#' @examples
#' textInput("caption", "Caption:", "Data Summary")
#' @export
textInput <- function(inputId, label, value = "") {
  tagList(
    label %AND% tags$label(label, `for` = inputId),
    tags$input(id = inputId, type="text", value=value)
  )
}
