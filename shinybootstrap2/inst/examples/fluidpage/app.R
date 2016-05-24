appExpr <- quote({
  shinyApp(
    ui = fluidPage(
      flowLayout(
        numericInput("rows", "How many rows?", 5),
        selectInput("letter", "Which letter?", LETTERS),
        sliderInput("value", "What value?", 0, 100, 50),
        selectInput("ui", "Input type", choices = c("numeric", "slider")),
        uiOutput("n_ui"),
        plotOutput("plot")
      )
    ),
    server = function(input, output) {
      output$n_ui <- renderUI({
        if (input$ui == "numeric")
          numericInput("n", "n", 1)
        else if (input$ui == "slider")
          sliderInput("n", "n", 1, 10, value = 1)
      })
      output$plot <- renderPlot( plot(head(cars, input$n)) )
    }
  )
})

# The useBS2 environment var controls whether or not we run the app using the
# wihtBootstrap2 wrapper function.
if (is.null(.GlobalEnv$useBS2))
  .GlobalEnv$useBS2 <- FALSE

if (.GlobalEnv$useBS2) {
  shinybootstrap2::withBootstrap2(appExpr, quoted = TRUE)
} else {
  eval(appExpr)
}
