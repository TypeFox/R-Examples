library(shiny)
library(catenary)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  cat <- reactive(function(){
    x <- runif(100,input$x0,input$x1)
    if(input$para_sim){
      y <- input$a * x^2 + input$b * x + input$c + rnorm(100,sd=input$sd)      
    } else {
      y <- f(x=x,c1=input$c1,c2=input$c2,lambda=input$lambda) + 
        rnorm(100,sd=input$sd)
    }
    cat <- fittedCatenary(x=x,y=y,R=input$R)
  })
  
  output$catPlot <- reactivePlot(function(){
    if(input$cat & input$para){
      p <- plot(cat(),fit=c("cat","para"),envelope=c("cat","para"))
    }
    if(!input$cat & input$para){
      p <- plot(cat(),fit=c("para"),envelope="para")
    }
    if(input$cat & !input$para){
      p <- plot(cat(),fit=c("cat"),envelope="cat")
    }
    if(!input$cat & !input$para){
      p <- plot(cat(),fit="none",envelope="none")
    }
    if(input$aspect){
      p <- p + coord_fixed()      
    }
    print(p)
  })
  
  output$summary <- reactivePrint(function(){
    Summary(cat())
  })

})