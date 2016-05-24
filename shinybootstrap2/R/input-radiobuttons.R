#' Create radio buttons
#'
#' Create a set of radio buttons used to select an item from a list.
#'
#' @param inputId Input variable to assign the control's value to
#' @param label Display label for the control, or \code{NULL}
#' @param choices List of values to select from (if elements of the list are
#' named then that name rather than the value is displayed to the user)
#' @param selected The initially selected value (if not specified then
#' defaults to the first value)
#' @param inline If \code{TRUE}, render the choices inline (i.e. horizontally)
#' @return A set of radio buttons that can be added to a UI definition.
#'
#' @family input elements
#' @seealso \code{\link{updateRadioButtons}}
#'
#' @examples
#' radioButtons("dist", "Distribution type:",
#'              c("Normal" = "norm",
#'                "Uniform" = "unif",
#'                "Log-normal" = "lnorm",
#'                "Exponential" = "exp"))
#' @export
radioButtons <- function(inputId, label, choices, selected = NULL, inline = FALSE) {
  # resolve names
  choices <- choicesWithNames(choices)

  # default value if it's not specified
  selected <- if (is.null(selected)) choices[[1]] else {
    validateSelected(selected, choices, inputId)
  }
  if (length(selected) > 1) stop("The 'selected' argument must be of length 1")

  options <- generateOptions(inputId, choices, selected, inline, type = 'radio')

  tags$div(id = inputId,
    class = 'control-group shiny-input-radiogroup',
    label %AND% tags$label(class = "control-label", `for` = inputId, label),
    options)
}
