library(shiny)
library(HH)

## Define server logic
shinyServer(function(input, output) {
  output$densityPlot <- renderPlot(
    {
      bivariateNormal(
        rho=input$rho,
        angle=input$angle,
        layout=c(1,1),
        colorkey=list(at=seq(0, .28, length=101)),
        zlim=c(0, .28)
      )
    })
})
