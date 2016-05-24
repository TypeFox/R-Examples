library(shiny)
library(catenary)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  cat <- reactive(function(){
    cat <- catenary(c1=input$c1,c2=input$c2,
                    lambda=input$lambda,x0=input$x0,x1=input$x1)
  })
  
  output$catPlot <- reactivePlot(function(){
    p <- plot(cat())
    print(p)
  })
  
  output$summary <- reactivePrint(function(){
    Summary(cat())
  })
  
})