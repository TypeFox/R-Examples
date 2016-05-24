# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)

shinyServer(function(input, output) {
  
  #output$distPlot <- renderPlot({
  
  # generate bins based on input$bins from ui.R
  # x <- faithful[, 2]
  #bins <- seq(min(x), max(x), length.out = input$bins + 1)
  
  # draw the histogram with the specified number of bins
  # hist(x, breaks = bins, col = 'darkgray', border = 'white')
  
  #})
  
  rng<-seq(-5,5,0.1)
  rngdf<-data.frame(rng)
  tplot<-ggplot(data=rngdf,aes(x=rng))+ylim(0,0.4)+
    ylab("Density")+xlab("Range of X values")+ggtitle("T density, with Normal(0,1) in black")
  
  rng2<-seq(-6,2,0.1)
  rng2df<-data.frame(rng2)
  tplot2<-ggplot(data=rng2df,aes(x=rng2))+ylim(0,0.07)+xlim(-6,-2)+
    ylab("Density")+xlab("Range of X values")+ggtitle("A closer view of the left tail")
  
  output$thePlot<-renderPlot({
    
    #make a normal distribution based on inputs.
    tplot+stat_function(fun=dt,args=list(df=input$degrees),color="steelblue")+stat_function(fun=dnorm,color="black")
  })
  
  output$secondPlot<-renderPlot({
    
    #make a normal distribution based on inputs.
    tplot2+stat_function(fun=dt,args=list(df=input$degrees),color="steelblue")+stat_function(fun=dnorm,color="black")
  })
  
})