#' Create a JSColor picker input widget
#' 
#' Creates a JSColor (Javascript/HTML Color Picker) widget to be used in shiny applications. 
#' 
#' @param inputId \code{\link{character}} (required): Specifies the input slot that will be used to access the value.
#' @param label \code{\link{character}}: Display label for the control, or NULL for no label.
#' @param value \code{\link{character}}: Initial RGB value of the color picker. Default is black ('#000000').
#' @param position \code{\link{character}}: Position of the picker relative to the text input ('bottom', 'left', 'top', 'right').
#' @param color \code{\link{character}}: Picker color scheme ('transparent' by default). Use RGB color coding ('000000').
#' @param mode \code{\link{character}}: Mode of hue, saturation and value. Can either be 'HSV' or 'HVS'.
#' @param slider \code{\link{logical}}: Show or hide the slider.
#' @param close \code{\link{logical}}: Show or hide a close button.
#'  
#' @seealso Other input.elements: \code{\link{animationOptions}}, \code{\link{sliderInput}}; 
#' \code{\link{checkboxGroupInput}}; \code{\link{checkboxInput}}; \code{\link{dateInput}}; 
#' \code{\link{dateRangeInput}}; \code{\link{fileInput}}; \code{\link{numericInput}}; 
#' \code{\link{passwordInput}}; \code{\link{radioButtons}}; \code{\link{selectInput}}, 
#' \code{\link{selectizeInput}}; \code{\link{submitButton}}; \code{\link{textInput}}
#'  
#' @examples 
#' # html code
#' jscolorInput("col", "Color", "21BF6B", slider = FALSE)
#' 
#' # example app
#' \dontrun{
#' shinyApp(
#' ui = fluidPage(
#'   jscolorInput(inputId = "col", label = "JSColor Picker", 
#'                value = "21BF6B", position = "right", 
#'                mode = "HVS", close = TRUE),
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
jscolorInput <- function(inputId, label, value, position = 'bottom', color = 'transparent', mode = 'HSV', slider = TRUE, close = FALSE) {
  tagList(        
    singleton(tags$head(tags$script(src = "RLumShiny/jscolor_inputBinding.js"))),
    singleton(tags$head(tags$script(src = "RLumShiny/jscolor/jscolor.js"))),
    if (missing(label)) { tags$p(" ") } else if (!is.null(label)) { tags$p(label) },
    tags$input(id = inputId, 
               value = ifelse(!missing(value), value, "#000000"), 
               class = sprintf("color {hash:true, pickerPosition:'%s', pickerBorderColor:'transparent', pickerFaceColor:'%s', pickerMode:'%s', slider:%s, pickerClosable:%s} shiny-bound-input", position, color, mode, tolower(slider), tolower(close)), 
               onchange = sprintf("$('#%s').trigger('afterChange')", inputId)),
    tags$script(sprintf("$('#%s').trigger('afterChange')", inputId))
  )
}