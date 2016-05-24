library(shiny)


#############################################
# Norma population info
############################################
muNorm <- 70
popMean <- muNorm
sigma <- 3

##############################
# Define server logic
##############################
shinyServer(function(input, output) {
  
  popDen <- reactive({
    xNorm <- seq(muNorm-5*sigma,muNorm+5*sigma,length.out=600)
    yNorm <- dnorm(xNorm,mean=muNorm,sd=sigma)
    list(x=xNorm,y=yNorm)
  })

  meanLineMax <- reactive({
    dnorm(popMean,mean=popMean,sd=sigma)
  })
  
  popMax <- reactive({
    den <- popDen()
    max(den$x)
  })
  
  popMin <- reactive({
    den <- popDen()
    min(den$x)
  })
  
  yMax <- reactive({
    densities <- popDen()
    max(densities$y)*1.5
  })
  
  sampleSize <- reactive({
    as.numeric(input$n)
  })
  
  
  output$go <- reactive({
    input$go
  })
  
  output$sizeSelect <- renderUI({
    if (input$actionType=="one") {
      selectInput(inputId="n",label="Sample Size",
                  choices=list("5"=5,"15"=15,"50"=50,"500"=500,"Ten Thousand"=10000),
                  selected="15")
    } else {
      sliderInput(inputId="n","Sample Size",min=2,max=50,value=10)
    }
  })

  
  output$actionType <- reactive({
    input$actionType
  })
  
  
  
  outputOptions(output, 'go', suspendWhenHidden=FALSE)
  outputOptions(output, 'actionType', suspendWhenHidden=FALSE)
  
 
  
  oneSample <- reactive({
    if (input$go > 0) {
    n <- isolate(sampleSize())
    sig <- sigma
    
    samp <- isolate(rnorm(n,mean=muNorm,sd=sig))
    samp
    } else c(69,71)
    
  })
  
  output$initialGraph <- renderPlot({
    popDen <- popDen()
    mu0 <- input$mu0
    ymax <- 1.2*meanLineMax()
    plot(popDen$x,popDen$y,type="l",lwd=3,col="red",
         main="Density Curve of Population",cex.main=1.5,
         xlab="",ylim=c(0,ymax),
         ylab="density")
    segments(x0=popMean,y0=0,x1=popMean,y1=meanLineMax(),
             lwd=4,col="red")
    if (mu0 != popMean) {
      segments(x0=mu0,y0=0,x1=mu0,y1=1.1*meanLineMax(),
               lwd=2,col="black")
      text(x=mu0,y=1.15*meanLineMax(),labels=expression(mu[0]),cex=1.5)
    }
    
  })
  
  output$initialGraph2 <- renderPlot({
    popDen <- popDen()
    mu0 <- input$mu0
    ymax <- 1.2*meanLineMax()
    plot(popDen$x,popDen$y,type="l",lwd=3,col="red",
         main="Density Curve of Population",cex.main=1.5,
         xlab="",ylim=c(0,ymax),
         ylab="density")
    segments(x0=popMean,y0=0,x1=popMean,y1=meanLineMax(),
             lwd=4,col="red")
    if (mu0 != popMean) {
      segments(x0=mu0,y0=0,x1=mu0,y1=1.1*meanLineMax(),
               lwd=2,col="black")
      text(x=mu0,y=1.15*meanLineMax(),labels=expression(mu[0]),cex=1.5)
    }
    
  })

  
  output$graphSample <- renderPlot({
    samp <- oneSample()
    popDen <- isolate(popDen())
    n <- isolate(sampleSize())
    mu0 <- isolate(input$mu0)
    
    xmin <- isolate(popMin())
    xmax <- isolate(popMax())
    ymax <- isolate(yMax())
    intLevel <- 0.80*ymax  #how high on plot to put conf interval
    
    conf = isolate(1-input$alpha)
    xbar = mean(samp)
    t.input = conf + ((1 - conf)/2)
    tMultiplier = qt(t.input, df = n - 1)
    se = sd(samp)/sqrt(n)
    margin = tMultiplier * se
    ci = c(xbar - margin, xbar + margin)
    
    hist(samp,freq=FALSE,col="lightblue",xlim=c(xmin,xmax),
         ylim=c(0,ymax),
         main="Population Density Curve\nand Histogram of Sample",
         xlab="",
         cex.main=2
    )
    lines(popDen,lwd=2,col="red")
    segments(x0=popMean,y0=0,x1=popMean,y1=isolate(meanLineMax()),
             lwd=4,col="red")
    segments(x0=mu0,y0=0,x1=mu0,y1=0.9*ymax,
             lwd=2,col="black")
    text(x=mu0,y=0.95*ymax,labels=expression(mu[0]),cex=1.5)
    
    
    segments(x0 = ci[1], y0 = intLevel, x1 = ci[2], y1 = intLevel, 
             col = "green", lwd = 3)
    text(x=ci[1],y=intLevel,labels="(")
    text(x=ci[2],y=intLevel,labels=")")
    points(xbar, intLevel, col = "blue", pch = 20,cex=2)
    
  })
  
  
  output$results <- renderTable({
    if (input$go > 0) {
      
    samp <- oneSample()
    n <- isolate(sampleSize())
    alpha <- isolate(input$alpha)
    xbar <- mean(samp)
    se <- sd(samp)/sqrt(n)
    
    mu0 <- isolate(input$mu0)
    tstat <- (xbar-mu0)/se
    
    pValue <- 2*pt(abs(tstat),df=n-1,lower.tail=FALSE)
    reject <- (pValue < alpha)
    
    nullCorrect <- (mu0 == popMean)
    
    error <- "none"
    
    if (nullCorrect && reject) error <- "Type-I"
    if ((!nullCorrect) && (!reject)) error <- "Type-II"
    
    results <- data.frame(
      NullHypothesis=paste("mu =",mu0),
      AltHypothesis=paste("mu !=",mu0),
      Pvalue=pValue,
      Decision=ifelse(reject,"Reject Null","Do Not Reject Null"),
      Error=error
      )
    results
    } #end if
    else data.frame(x=1)  #just something to return before display is possible
    
  })
  
  intervalFrame <- reactive({
    
    input$go
    n <- isolate(input$n)
    mu0 <- isolate(input$mu0)
    alpha <- isolate(input$alpha)
    action <- isolate(input$actionType)
    
    if (action == "fiveThousand") {
      itemNumb <- 5000*n
      sampleItems <- rnorm(itemNumb,mean=muNorm,sd=sigma)
      
      sampleMatrix <- matrix(sampleItems,ncol=n,nrow=5000)
      
      xbar <- rowSums(sampleMatrix)/n
      se <- sqrt((rowSums(sampleMatrix^2)-n*xbar^2)/(n^2-n))
      
      tstat <- (xbar-mu0)/se
      
      pValue <- 2*pt(abs(tstat),df=n-1,lower.tail=FALSE)
      reject <- (pValue < alpha)
      
      nullCorrect <- (mu0 == popMean)
      error <- character(5000)
      
      if (nullCorrect) {
        error <- ifelse(reject,"Type-I Error","none")
      } else {
        error <- ifelse(reject,"none","Type-II Error")
      }
    
      results <- data.frame(
        SampleMean=xbar,
        tStatistic=tstat,
        TwoSidedPvalue=pValue,
        Decision=ifelse(reject,"Reject Null","Do Not Reject Null"),
        Error=error
      )
      results
    }
  })
  
  output$summary <- renderTable({
    
    frame <- intervalFrame()
    
    nullCorrect <- isolate(input$mu0 == popMean)
    OKCount <- length(frame$Error[frame$Error=="none"])
    
    if (nullCorrect) {
      tab <- data.frame("Correct Decisions"=as.character(OKCount),
                      "Type-I Errors"=as.character(5000-OKCount),
                      "Type-I Error Rate"=paste(round(100-OKCount/50,2),"%",sep=""))
      tab
    } else {
      tab <- data.frame("Correct Decisions"=as.character(OKCount),
                        "Type-II Errors"=as.character(5000-OKCount),
                        "Type-II Error Rate"=paste(round(100-OKCount/50,2),"%",sep=""))
      tab
    }
    
  })
  
  output$intervalFrame = renderDataTable({
    intervalFrame()[,2:5]}, options = list(aLengthMenu = c(5, 30, 50), iDisplayLength = 5)
  )

  
})
  