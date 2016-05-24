library(shiny)
library(lattice)

## Define server logic
shinyServer(function(input, output) {

  dataset <- reactive({
    set.seed(input$seed)
    data.frame(x=rnorm(100),
               e=rnorm(100))
  })

  y <- reactive({
    dataset()$x*input$rho + dataset()$e*(1-input$rho^2)^.5
  })

  output$correlationPlot <- renderPlot(
    {
      maxabs <- c(-1,1) * 3.9
      xyplot(y() ~ dataset()$x, aspect="iso",
             xlim=maxabs, ylim=maxabs, scales=list(at=c(-2,0,2)),
             xlab="x", ylab=list("y", rot=0),
             main=as.expression(substitute(rho == r, c(alist(rho=rho), list(r=input$rho)))),
             pch=19, cex=1.3, col="blue")
    })
})
