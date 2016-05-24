#' Change the value of a slider input on the client
#'
#' @template update-input
#' @param value The value to set for the input object.
#'
#' @seealso \code{\link{sliderInput}}
#'
#' @examples
#' \dontrun{
#' shinyServer(function(input, output, session) {
#'
#'   observe({
#'     # We'll use the input$controller variable multiple times, so save it as x
#'     # for convenience.
#'     x <- input$controller
#'
#'     # Similar to number and text. only label and value can be set for slider
#'     updateSliderInput(session, "inSlider",
#'       label = paste("Slider label", x),
#'       value = x)
#'
#'     # For sliders that pick out a range, pass in a vector of 2 values.
#'     updateSliderInput(session, "inSlider2", value = c(x-1, x+1))
#'
#'     # An NA means to not change that value (the low or high one)
#'     updateSliderInput(session, "inSlider3", value = c(NA, x+2))
#'   })
#' })
#' }
#' @export
updateSliderInput <- function(session, inputId, label = NULL, value = NULL) {
  message <- dropNulls(list(label=label, value=value))
  session$sendInputMessage(inputId, message)
}

updateInputOptions <- function(session, inputId, label = NULL, choices = NULL,
                               selected = NULL, inline = FALSE,
                               type = 'checkbox') {

  choices <- choicesWithNames(choices)
  if (!is.null(selected))
    selected <- validateSelected(selected, choices, inputId)

  options <- if (length(choices))
    format(tagList(
      generateOptions(inputId, choices, selected, inline, type = type)
    ))

  message <- dropNulls(list(label = label, options = options, value = selected))

  session$sendInputMessage(inputId, message)
}

#' Change the value of a checkbox group input on the client
#'
#' @template update-input
#' @inheritParams checkboxGroupInput
#'
#' @seealso \code{\link{checkboxGroupInput}}
#'
#' @examples
#' \dontrun{
#' shinyServer(function(input, output, session) {
#'
#'   observe({
#'     # We'll use the input$controller variable multiple times, so save it as x
#'     # for convenience.
#'     x <- input$controller
#'
#'     # Create a list of new options, where the name of the items is something
#'     # like 'option label x 1', and the values are 'option-x-1'.
#'     cb_options <- list()
#'     cb_options[[sprintf("option label %d 1", x)]] <- sprintf("option-%d-1", x)
#'     cb_options[[sprintf("option label %d 2", x)]] <- sprintf("option-%d-2", x)
#'
#'     # Change values for input$inCheckboxGroup
#'     updateCheckboxGroupInput(session, "inCheckboxGroup", choices = cb_options)
#'
#'     # Can also set the label and select items
#'     updateCheckboxGroupInput(session, "inCheckboxGroup2",
#'       label = paste("checkboxgroup label", x),
#'       choices = cb_options,
#'       selected = sprintf("option-%d-2", x)
#'     )
#'   })
#' })
#' }
#' @export
updateCheckboxGroupInput <- function(session, inputId, label = NULL,
                                     choices = NULL, selected = NULL,
                                     inline = FALSE) {
  updateInputOptions(session, inputId, label, choices, selected, inline)
}


#' Change the value of a radio input on the client
#'
#' @template update-input
#' @inheritParams radioButtons
#'
#' @seealso \code{\link{radioButtons}}
#'
#' @examples
#' \dontrun{
#' shinyServer(function(input, output, session) {
#'
#'   observe({
#'     # We'll use the input$controller variable multiple times, so save it as x
#'     # for convenience.
#'     x <- input$controller
#'
#'     r_options <- list()
#'     r_options[[sprintf("option label %d 1", x)]] <- sprintf("option-%d-1", x)
#'     r_options[[sprintf("option label %d 2", x)]] <- sprintf("option-%d-2", x)
#'
#'     # Change values for input$inRadio
#'     updateRadioButtons(session, "inRadio", choices = r_options)
#'
#'     # Can also set the label and select an item
#'     updateRadioButtons(session, "inRadio2",
#'       label = paste("Radio label", x),
#'       choices = r_options,
#'       selected = sprintf("option-%d-2", x)
#'     )
#'   })
#' })
#' }
#' @export
updateRadioButtons <- function(session, inputId, label = NULL, choices = NULL,
                               selected = NULL, inline = FALSE) {
  # you must select at least one radio button
  if (is.null(selected) && !is.null(choices)) selected <- choices[[1]]
  updateInputOptions(session, inputId, label, choices, selected, inline, type = 'radio')
}
