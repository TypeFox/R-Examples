library(shiny)
library(catenary)


# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  cat <- reactive(function(){
    endpoints <- matrix(c(input$x0,input$x1,input$y0,input$y1),ncol=2)
    if(input$natural=="nat"){
      cat <- catenary(endpoints=endpoints,type='natural')
    } else if( input$natural == "max"){
      cat <- catenary(endpoints=endpoints,type='max')
    } else {
      cat <- catenary(endpoints=endpoints,L = input$L)
    }
    return(cat)
  })
  
  output$catPlot <- reactivePlot(function(){
    p <- plot(cat())
    print(p)
  })
  
  output$lengthslider <- reactiveUI(function() {
    endpoints <- matrix(c(input$x0,input$x1,input$y0,input$y1),ncol=2)
    tmp <- minmaxLength(endpoints=endpoints)
    xmin <- tmp[1]
    xmax <- tmp[2]
    
    sliderInput(inputId = "L",
                label = paste("Limit range"),
                min = xmin, max = xmax, value = mean(c(xmin,xmax)))
  })
  

})