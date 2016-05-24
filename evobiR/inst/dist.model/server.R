library(shiny)


shinyServer(function(input, output) {
  
  x <- reactive({
    if(input$select == 1){
      seq(input$mu - 40, 
          input$mu + 40, 
          length = 200)
    }else if(input$select == 2){
      high <- qexp(.999, rate = input$lambda, lower.tail = TRUE, log.p = FALSE)
      seq(0, high, length = 200)
    }else if(input$select == 3){
      high <- qgamma(.999, shape = input$kappa, rate = input$theta)
      seq(0, high, length = 200)
    }else if(input$select == 4){
      seq(input$mu2 - 30, 
          input$mu2 + 30, 
          length = 200)
    }else if(input$select == 5){
      seq(1:50)
    }else if(input$select == 6){
      seq(0, 1, length = 200)
    }
  })
  y <- reactive({
    if(input$select == 1){
      dnorm(x(), mean = input$mu, sd=input$sigma)
    }else if(input$select == 2){
      dexp(x(), rate = input$lambda)
    }else if(input$select == 3){
      dgamma(x(), shape = input$kappa, rate = input$theta)
    }else if(input$select == 4){
      dlogis(x(), location = input$mu2, scale = input$sigma2*0.551328895)
    }else if(input$select == 5){
      dpois(x(), lambda = input$lambda2)
    }else if(input$select == 6){
      dbeta(x(), shape1 = input$alpha, shape2= input$beta)
    }
})

SE <- reactive({
  if(input$select == 1){
    sd(rnorm(input$n, mean = input$mu, sd=input$sigma)) / sqrt(input$n)
  }
})
SD <- reactive({
  if(input$select == 1){
    sd(rnorm(input$n, mean = input$mu, sd=input$sigma))
  }
})



  output$treePlot <- renderPlot({
    if(input$select == 1){
      plot(x=x(), y=y(), col = "red", ylab="density", ylim=c(0,.6), xlab="x",
           main=paste("Normal probability density"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
      mtext(text = paste("SE =", round(SE(), digits=4)),
            side=3,line=-2)
      mtext(text = paste("SD est =", round(SD(), digits=4)),
            side=3,line=-3)
    }else if(input$select == 2){
      plot(x=x(), y=y(), col = "red", ylab="density", ylim=c(0,11), xlab="x",
           main=paste("Exponential probability density"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
    }else if(input$select == 3){
      plot(x=x(), y=y(), col = "red", ylab="density",  xlab="x",
           main=paste("Gamma probability density"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
    }else if(input$select == 4){
      plot(x=x(), y=y(), col = "red", ylab="density",  xlab="x",
           main=paste("Logistic probability density"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
    }else if(input$select == 5){
      plot(x=x(), y=y(), col = "red", ylab="density",  xlab="x",
           main=paste("Poisson probability density"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
      points(x=x(), y=y(), col="blue", pch=19)
    }else if(input$select == 6){
      plot(x=x(), y=y(), col = "red", ylab="density",  xlab="x",
           main=paste("Beta probability density"), type="l", lwd=3, ylim=c(0,6))
      abline(h=0, lty=3, cex=2)
    }
  })  
})

