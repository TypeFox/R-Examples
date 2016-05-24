
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
   # x    <- faithful[, 2]
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
  #  hist(x, breaks = bins, col = 'darkgray', border = 'white')

  #})
  
  
  output$binomialPlot<-renderPlot({
    
    #make a binomial density based on inputs.
    rng=seq(0,input$bothn,1)
    values<-dbinom(x=rng,p=input$bothp,size=input$bothn)
    bothdf<-data.frame(rng,values)
    
    ggplot(bothdf)+geom_bar(aes(y=values,x=rng),stat="identity",fill="navy blue",binwidth=1/length(rng))+
      ylab("Probability of X heads")+xlab("Number of Heads")+
      ggtitle("Binomial Probability Plot") +coord_cartesian(ylim=c(0,1),xlim=c(-0.5,50.5))
  })
  
  output$pbinomialPlot<-renderPlot({
    
    #make a binomial density based on inputs.
    rng=seq(0,10,1)
    values<-dbinom(x=rng,p=input$singlep,size=10)
    bothdf<-data.frame(rng,values)
    
    ggplot(bothdf)+geom_bar(aes(y=values,x=rng),stat="identity",fill="navy blue",binwidth=1/length(rng))+
      ylab("Probability of X heads")+xlab("Number of Heads")+
      ggtitle("Binomial Probability Plot") +coord_cartesian(ylim=c(0,1),xlim=c(-0.5,10.5))
  })
  
  output$singletext <- renderText({ 
    paste("Our expected number of heads is 10 *",input$singlep,", which is ", 10*input$singlep,". We don't round this to the nearest value, by the way.")
  })
  
  output$bothtext <- renderText({ 
    paste("Our expected number of heads is n * p, or ",input$bothn,"*",input$bothp,"=",input$bothn*input$bothp,". We don't round this to the nearest value, by the way.")
  })
  
})
