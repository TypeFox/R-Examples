library(shiny)
library(catenary)

data(ctesiphon)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  cat <- reactive(function(){
    if(input$internal=="int"){
      cat <- fittedCatenary(x=ctesiphon$internal$x,y=ctesiphon$internal$y,R=10)
    } else {
      cat <- fittedCatenary(x=ctesiphon$external$x,y=ctesiphon$external$y,R=10)
    }
    return(cat)
  })
  
  output$catPlot <- reactivePlot(function(){
    p <- plot(cat())
    print(p)
  })
  
  output$summary <- reactivePrint(function(){
    Summary(cat())
  })
})