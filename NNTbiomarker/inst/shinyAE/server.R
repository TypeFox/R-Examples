appName = "shinyAE"
cat("Launching ", appName, "\n")

library(shiny)

# options(shiny.error=recover)

shinyServer(function(input, output, session) {

  thisSession <<- session

  source("skislope.R", local=TRUE)
  source("AEplot.R", local=TRUE)

  rV = reactiveValues(RSclicked = 30)

  observe({
    input$ODXlow
    updateNumericInput(inputId = "RSchosen",
                       session = thisSession,
                       value = OncotypeRScutoffs[1])
  })
  observe({
    input$ODXhigh
    updateNumericInput(inputId = "RSchosen",
                       session = thisSession,
                       value = OncotypeRScutoffs[2])
  })
  observe({
    input$TailorX_low
    updateNumericInput(inputId = "RSchosen",
                       session = thisSession,
                       value = TailorXRScutoffs[1])
  })
  observe({
    input$TailorX_high
    updateNumericInput(inputId = "RSchosen",
                       session = thisSession,
                       value = TailorXRScutoffs[2])
  })

  output$skislope = renderPlot({
#     if(is.null(input$skislope_click))
#       RSarg = 49
#     else RSarg = input$skislope_click$x
#     RSarg = round(max(RSarg, 5))
#     cat("============  Setting rV$RSclicked: ", rV$RSclicked,
#         " ======\n")
#     isolate({
#       rV$RSclicked = RSarg
#     })
    skisloplot(input$RSchosen, input$ytop)
  })

  output$AEstack = renderPlot({
    AEplot(input$RSchosen, makeTitle=TRUE)
  })

#   observe({
#     RSclicked = try(input$skislope_click$x)
#     if(class(RSclicked) == "try-error")
#       rV$RSclicked = 50
#     else
#       rV$RSclicked = RSclicked
#   })
})
