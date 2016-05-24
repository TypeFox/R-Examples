
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)


perm = function(n, x) {
  return(factorial(n) / factorial(n-x))
}

comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}

binlik = function(n, s, p) {
  return (comb(n,s)*p^(s)*(1-p)^(n-s))
}

ourcolor="#0A5E4F"

shinyServer(function(input, output) {

  
  
  output$betapriorPlot <- renderPlot({
    x=seq(0,1,0.01)
    y=dbeta(x=x, shape1 = input$alpha, shape2 = input$beta)
    betapriordf <- data.frame(cbind(x,y))
    ggplot(betapriordf) + geom_ribbon(aes(x = x, ymax = y,ymin=0),fill=ourcolor,color=ourcolor)+
      xlab("Value of p")+ylab("Density")+
      ggtitle(paste("Beta(",input$alpha,",",input$beta,") Prior",sep=""))
    
  })
  
  
  
  output$binlikPlot <- renderPlot({
  
    x = seq(0, 1, 0.01)
    y = binlik(input$n, input$s, x)
    likdf <- data.frame(cbind(x, y))
    ggplot(likdf) + geom_ribbon(aes(x = x, ymax = y,ymin=0),fill=ourcolor,color=ourcolor)+
      xlab("Value of p")+ylab("Likelihood of that value")+
      ggtitle(paste("Likelihood of p for Bin(",input$n,",",round(input$s/input$n,2),")",sep=""))
  
  
  })
  
  

  
  
  output$betapostPlot <- renderPlot ({
    x=seq(0,1,0.01)
    y=dbeta(x=x, shape1 = input$alpha + input$s, shape2 = input$beta + input$n - input$s)
    betapriordf <- data.frame(cbind(x,y))
    ggplot(betapriordf) + geom_ribbon(aes(x = x, ymax=y,ymin=0),fill=ourcolor,color=ourcolor)+
      xlab("Value of p")+ylab("Density")+
      ggtitle(paste("Beta (",input$alpha + input$s,",",input$beta + input$n - input$s,") Posterior",sep=""))
  })
  
  
  
  output$sUI <- renderUI({ 
    sliderInput("s",
                "Number of successes for the binomial:",
                min = 0,
                max = as.numeric(input$n),
                value = 1,
                step=1)
  })
  
  
})
