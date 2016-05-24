#' Create a bootstrap tooltip
#' 
#' Create bootstrap tooltips for any HTML element to be used in shiny applications. 
#' 
#' @param refId \code{\link{character}} (required): id of the element the tooltip is to be attached to.
#' @param text \code{\link{character}}: Text to be displayed in the tooltip.
#' @param attr \code{\link{character}}: Attach tooltip to all elements with attribute \code{attr='refId'}.
#' @param animation \code{\link{logical}}: Apply a CSS fade transition to the tooltip.
#' @param delay \code{\link{numeric}}: Delay showing and hiding the tooltip (ms).
#' @param html \code{\link{logical}}: Insert HTML into the tooltip.
#' @param placement \code{\link{character}}: How to position the tooltip - top | bottom | left | right | auto. When 'auto' is specified, it will dynamically reorient the tooltip. For example, if placement is 'auto left', the tooltip will display to the left when possible, otherwise it will display right.
#' @param trigger \code{\link{character}}: How tooltip is triggered - click | hover | focus | manual. You may pass multiple triggers; separate them with a space.
#' 
#' @examples 
#' # javascript code
#' tt <- tooltip("elementId", "This is a tooltip.")
#' str(tt)
#' 
#' # example app
#' \dontrun{
#' shinyApp(
#' ui = fluidPage(
#'   jscolorInput(inputId = "col", label = "JSColor Picker", 
#'                value = "21BF6B", position = "right", 
#'                mode = "HVS", close = TRUE),
#'   tooltip("col", "This is a JScolor widget"),
#'   
#'   checkboxInput("cbox", "Checkbox", FALSE),
#'   tooltip("cbox", "This is a checkbox"),
#'   
#'   checkboxGroupInput("cboxg", "Checkbox group", selected = "a", 
#'                      choices = c("a" = "a",
#'                                  "b" = "b",
#'                                  "c" = "c")),
#'   tooltip("cboxg", "This is a <b>checkbox group</b>", html = TRUE),
#'   
#'   selectInput("select", "Selectinput", selected = "a", choices = c("a"="a", "b"="b")),
#'   tooltip("select", "This is a text input field", attr = "for", placement = "right"),
#'   
#'   passwordInput("pwIn", "Passwordinput"),
#'   tooltip("pwIn", "This is a password input field"),
#'   
#'   plotOutput("plot")
#' ),
#' server = function(input, output) {
#'   output$plot <- renderPlot({
#'     plot(cars, col = input$col, cex = 2, pch = 16)
#'  })
#' })
#' }
#' @import shiny
#' @export
tooltip <- function(refId, text, attr = NULL, animation = TRUE, delay = 100, html = TRUE, placement = 'auto', trigger = 'hover') {
  if (is.null(attr))
    el <- sprintf("'#%s'", refId)
  else
    el <- sprintf("\"[%s='%s']\"", attr, refId)
  
  tagList(
      tags$head(tags$script(HTML(sprintf("$(window).load(function(){ $(%s).tooltip({ html: %s, trigger: '%s', title: '%s', animation: %s, delay: {'show': %i, 'hide': %i}, placement: '%s' }); })",
                                  el, tolower(html), trigger, text, tolower(animation), delay, delay, placement))))
  )
}