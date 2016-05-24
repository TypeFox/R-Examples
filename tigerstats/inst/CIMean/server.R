library(shiny)
library(scales) # for transparency in density-plot fill

## Set up underlying populations:
source("setup.R")

## Set upper limit on sims
simLimit <- 10000 #upper limit on number of sims at once


########################################
## the server
#########################################

function(input, output, session) {
  ## set see so that users arelikely to get different results
  set.seed(as.numeric(Sys.time()))
  
  ########################################## 
  ## server material for coverage tab:
  ##########################################
  
  rv <- reactiveValues(
    popDen = normalDen,
    popMean = normalMean,
    popMax = max(normalDen$x),
    popMin = min(normalDen$x),
    yMax = 1.5*max(normalDen$y),
    sample = NULL, 
    mean = NULL, 
    lower = NULL,
    upper = NULL,
    sims = 0,
    good = 0,
    begin = TRUE,
    tstats = numeric())
  
  # respond to request to change population
  observeEvent(
    input$popDist,
    {
    rv$popDen <- switch(input$popDist,
                        normal=normalDen,
                        skew=skewDen,
                        superskew=superSkewDen,
                        outliers=outlierDen)
    rv$popMean <- switch(input$popDist,
                         normal=normalMean,
                         skew=skewMean,
                         superskew=superSkewMean,
                         outliers=outlierMean)
    rv$popMax <- switch(input$popDist,
                        normal=max(normalDen$x),
                        skew=max(skewDen$x),
                        superskew=max(superSkewDen$x),
                        outliers=max(outlierDen$x))
    rv$popMin <- switch(input$popDist,
                        normal=min(normalDen$x),
                        skew=min(skewDen$x),
                        superskew=min(superSkewDen$x),
                        outliers=min(outlierDen$x))
    rv$yMax <- switch(input$popDist,
                      normal=1.5*max(normalDen$y),
                      skew=1.5*max(skewDen$y),
                      superskew=1.5*max(superSkewDen$y),
                      outliers=1.5*max(outlierDen$y))
    }
  )
  
  # respond to request for sample(s)
  observeEvent(
    input$takeSample, 
    {
    # if user is at the beginning or has started over, and one sample at a time is
    # requested, then we will assume that the user would prefer to stat on the "Latest Interval"
    # tabPanel, rather than the "t-statistic" tabPanel.  Assure this:
    if (rv$begin && input$sims == 1) {
      updateTabsetPanel(session, "coverageTabsetPanel", selected = "Latest Interval")
      }
    # get the samples, make the intervals
    n <- input$n
    reps <- min(input$sims, simLimit)
    # grab all the random items you need at once:
    itemNumb <- reps*n
    sampleItems <- switch(input$popDist,
                          normal=rnorm(itemNumb,mean=muNorm,sd=sigmaNorm),
                          skew=rgamma(itemNumb,shape=shapeGamma,scale=scaleGamma),
                          superskew=rpareto(itemNumb,alpha=alphaPareto,theta=thetaPareto),
                          outliers=routlier(itemNumb)
                          )
    # arrange the random items in a matrix; the rows are your samples
    sampleMatrix <- matrix(sampleItems,ncol=n,nrow=reps)
    conf = input$confLevel/100
    t.input = conf + ((1 - conf)/2)
    tMultiplier = qt(t.input, df = n - 1)
    # from the matrix, quickly compute the items you need
    xbar <- rowSums(sampleMatrix)/n
    se <- sqrt((rowSums(sampleMatrix^2)-n*xbar^2)/(n^2-n))
    margin = tMultiplier * se
    lower <- xbar - margin
    upper <- xbar + margin
    goodInterval <- ((rv$popMean > lower) & (rv$popMean < upper))
    goodCount <- sum(goodInterval)
    latestSamp <<- sampleMatrix[reps,]
    # store in rv
    rv$sample <- sampleMatrix[reps, ]
    rv$mean <- xbar[reps]
    rv$lower <- lower[reps]
    rv$upper <- upper[reps]
    rv$sims <- rv$sims + reps
    rv$good <- rv$good + goodCount
    rv$begin <- FALSE
    rv$tstats <- c(rv$tstats, (xbar-rv$popMean)/se)
    }
    )
  
  # respond to request to start over
  observeEvent(input$reset,
               {
                 rv$sample <- NULL
                 rv$mean <- NULL
                 rv$lower <- NULL
                 rv$upper <- NULL
                 rv$sims <- 0
                 rv$good <- 0
                 rv$begin <- TRUE
                 rv$tstats <- numeric()
               })
  
  # need output to help with conditional panels
  output$beginning <- reactive({
    rv$begin
  })
  
  # needed for the conditional panels to work
  outputOptions(output, 'beginning', suspendWhenHidden=FALSE)
  
  # now for output that the user will see
  output$initialGraph <- renderPlot({
    # the underlying population
    plot(rv$popDen$x,rv$popDen$y,type="l",lwd=3,col="red",
         main="Density Curve of Population",
         xlim=c(rv$popMin,rv$popMax),
         ylim=c(0,rv$yMax),
         xlab="",
         ylab="density")
    abline(v=rv$popMean,lwd=2)
  })
  
  output$plotSample <- renderPlot({
    # the underlying population
    plot(rv$popDen$x,rv$popDen$y,type="l",lwd=3,col="red",
         main="Density Curve of Population, with Random Sample",
         xlim=c(rv$popMin,rv$popMax),
         ylim=c(0,rv$yMax),
         xlab="",
         ylab="density")
    abline(v=rv$popMean,lwd=2)
    # sample and interval
    if (! rv$begin) {
      # density plot for the sample
      sampDen <- density(rv$sample, from = 0)
      xdens <- sampDen$x
      ydens <- sampDen$y
      firstx <- xdens[1]
      lastx <- xdens[length(xdens)]
      polygon(x = c(firstx,xdens,lastx), y = c(0,ydens,0), col = alpha("lightblue",0.5))
      # now the interval
      intLevel <- 0.95*rv$yMax
      segments(x0 = rv$lower, y0 = intLevel, x1 = rv$upper, y1 = intLevel, 
               col = "green", lwd = 3)
      text(x=rv$lower,y=intLevel,labels="(")
      text(x=rv$upper,y=intLevel,labels=")")
      points(rv$mean, intLevel, col = "blue", pch = 20,cex=2)
      rug(rv$sample)
      }
    })
  # summary of intervals so far
  output$summary <- renderTable({
    df <- data.frame(rv$sims,
                     rv$good,
                     ifelse(rv$sims >0, round(rv$good/rv$sims*100,3), NA))
    names(df) <- c("Simulations", "Good Intervals", "Percentage Good")
    df
    },
    include.rownames = FALSE)
  # produce the t-statistic plot
  output$tstatistic <- renderPlot({
    input$takeSample
    n <- input$n
    numberSims <- rv$sims
    tstats <- rv$tstats
    if (numberSims == 1) {
      tstatDen <- density(tstats,n=1024,from=-10,to=10,bw=1)
      }
    if (numberSims >= 2 && n < 5) {
      tstatDen <- density(tstats,n=1024,from=-10,to=10,bw=0.1)
      }
    if (numberSims >= 2 && n >= 5) {
      tstatDen <- density(tstats,n=1024,from=-10,to=10,bw="SJ")
      }
    if (numberSims > 0) {
      ymax <- max(tstatDen$y,dt(0,df=n-1))
      plot(tstatDen$x,tstatDen$y,type="l",lwd=2,col="blue",
           main="t-statistic vs. t-curve",cex.main=2,
           xlab="t", ylim=c(0,ymax),xlim=c(-6,6),
           ylab="density")
      curve(dt(x,df=n-1),-6,6,col="red",lwd=2,add=TRUE)
      } #end check that there are samples
  })
  
  #####################################################  
  ## server material for fifty-at-a-time tab
  #####################################################
  
  rv2 <- reactiveValues(lower = NULL,
                        upper = NULL,
                        good = NULL,
                        low = NULL,
                        high = NULL,
                        begin = TRUE,
                        popDen = normalDen,
                        popMean = normalMean,
                        popMax = max(normalDen$x),
                        popMin = min(normalDen$x),
                        yMax = 1.5*max(normalDen$y))
  
  observeEvent(
    input$popDist2,
    {
    rv2$popDen <- switch(input$popDist2,
                         normal=normalDen,
                         skew=skewDen,
                         superskew=superSkewDen,
                         outliers=outlierDen)
    rv2$popMean <- switch(input$popDist2,
                          normal=normalMean,
                          skew=skewMean,
                          superskew=superSkewMean,
                          outliers=outlierMean)
    rv2$popMax <- switch(input$popDist2,
                         normal=max(normalDen$x),
                         skew=max(skewDen$x),
                         superskew=max(superSkewDen$x),
                         outliers=max(outlierDen$x))
    rv2$popMin <- switch(input$popDist2,
                         normal=min(normalDen$x),
                         skew=min(skewDen$x),
                         superskew=min(superSkewDen$x),
                         outliers=min(outlierDen$x))
   rv2$yMax <- switch(input$popDist2,
                      normal = 1.5*max(normalDen$y),
                      skew = 1.5*max(skewDen$y),
                      superskew = 1.5*max(superSkewDen$y),
                      outliers = 1.5*max(outlierDen$y))
   }
   )

  observeEvent(
    input$takeSample2,
    {
    # get the samples, make the intervals
    n <- input$n2
    reps <- 50
    # grab all the random items you need at once:
    itemNumb <- reps*n
    sampleItems <- switch(input$popDist2,
                          normal=rnorm(itemNumb,mean=muNorm,sd=sigmaNorm),
                          skew=rgamma(itemNumb,shape=shapeGamma,scale=scaleGamma),
                          superskew=rpareto(itemNumb,alpha=alphaPareto,theta=thetaPareto),
                          outliers=routlier(itemNumb))
    # arrange the random items in a matrix; the rows are your samples
    sampleMatrix <- matrix(sampleItems,ncol=n,nrow=reps)
    conf = input$confLevel2/100
    t.input = conf + ((1 - conf)/2)
    tMultiplier = qt(t.input, df = n - 1)
    # from the matrix, quickly compute the items you need
    xbar <- rowSums(sampleMatrix)/n
    se <- sqrt((rowSums(sampleMatrix^2)-n*xbar^2)/(n^2-n))
    margin = tMultiplier * se
    lower <- xbar - margin
    upper <- xbar + margin
    highInterval <- (rv2$popMean < lower)
    lowInterval <- (rv2$popMean > upper)
    goodInterval <- !(lowInterval | highInterval)
    # store in rv2
    rv2$lower <- lower
    rv2$upper <- upper
    rv2$good <- goodInterval
    rv2$low <- lowInterval
    rv2$high <- highInterval
    rv2$begin <- FALSE
    }
    )
  
  observeEvent(
    input$reset2,
    {
    rv2$lower <- NULL
    rv2$upper <- NULL
    rv2$good <- NULL
    rv2$low <- NULL
    rv2$high <- NULL
    rv2$begin <- TRUE
    }
    )
  
  output$beginning2 <- reactive({
    rv2$begin
  })
  
  # needed for the conditional panels to work
  outputOptions(output, 'beginning2', suspendWhenHidden=FALSE)
  
  output$initialGraph2 <- renderPlot({
    # the underlying population
    plot(rv2$popDen$x,rv2$popDen$y,type="l",lwd=3,col="red",
         main="Density Curve of Population",
         xlim=c(rv2$popMin,rv2$popMax),
         ylim=c(0,rv2$yMax),
         xlab="",
         ylab="density")
    abline(v=rv2$popMean,lwd=2)
  })
  
  output$plotSample2 <- renderPlot({
    # the underlying population
    plot(rv2$popDen$x,rv2$popDen$y,type="l",lwd=3,col="red",
         main="Density Curve of Population, with Intervals",
         xlim=c(rv2$popMin,rv2$popMax),
         ylim=c(0,rv2$yMax),
         xlab="",
         ylab="density")
    abline(v=rv2$popMean,lwd=2)
    # intervals
    if (! rv2$begin) {
      reps <- length(rv2$lower)
      for(i in 1:reps) {
        interval <- c(rv2$lower[i],rv2$upper[i])
        color <- ifelse(rv2$good[i],"green","red")
        width <- ifelse(rv2$good[i],1,2)
        height <- 0.98*rv2$yMax *i/50
        lines(interval, c(height,height), col = color, lwd = width)
      }
      }
    })
  
  # summary of intervals so far
  # argument include.rownames is passed to xtable
  output$summary2 <- renderTable({
    df <- data.frame(50,
                     sum(rv2$good),
                     sum(rv2$low),
                     sum(rv2$high))
    names(df) <- c("Simulations", "Good Intervals", "Low Intervals", "High Intervals")
    df}, include.rownames = FALSE)
  } # at long last, we end the server function