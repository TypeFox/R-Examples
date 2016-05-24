library(shiny)

########################################
#  Utility Functions for Generating Pareto Values
########################################

rpareto <- function(n,alpha,theta) {#random values for Pareto(alpha,theta) distribution
  theta*((1-runif(n))^(-1/alpha)-1)
}

dpareto <- function(x,alpha,theta) {  #pdf for Pareto(alpha,theta) distribution
  alpha*theta^alpha/(x+theta)^(alpha+1)
}


#############################################
# Generate the populations
############################################
muNorm <- 70
sigmaNorm <- 5
shapeGamma <- 2
scaleGamma <- 50

# for pareto:
alphaPareto <- 5
thetaPareto <- 100
tailProb <- 0.02  #want to find a Value at risk of 1 - this
valRisk <- thetaPareto*(tailProb^(-.5)-1)

# for pop with group of outliers
propOutliers <- 0.10
meanOutliers <- 200
sdOutliers <- 5
meanRegulars <- 50
sdRegulars <- 5

routlier <- function(n) {
  propNormals <- 1- propOutliers
  whichHump <- rbinom(n,size=1,prob=propNormals)
  outlierSamp <- ifelse(whichHump,rnorm(n,mean=meanRegulars,sd=sdRegulars),
                        rnorm(n,mean=meanOutliers,sd=sdOutliers))
  outlierSamp
}

######################################
# Make population densities
#####################################
xNorm <- seq(muNorm-5*sigmaNorm,muNorm+5*sigmaNorm,length.out=600)
yNorm <- dnorm(xNorm,mean=muNorm,sd=sigmaNorm)
normalDen <- list(x=xNorm,y=yNorm)

xSkew <- seq(0,shapeGamma*scaleGamma+7.5*sqrt(shapeGamma)*scaleGamma,
             length.out=600)
ySkew <- dgamma(xSkew,shape=shapeGamma,scale=scaleGamma)
skewDen <- list(x=xSkew,y=ySkew)

xSuperSkew <- seq(0,valRisk,length.out=600)
ySuperSkew <- dpareto(xSuperSkew,alpha=alphaPareto,theta=thetaPareto)
superSkewDen <- list(x=xSuperSkew,y=ySuperSkew)

xOut <- seq(0,meanOutliers+5*sdOutliers,length.out=600)
yOut <- (1-propOutliers)*dnorm(xOut,mean=meanRegulars,sd=sdRegulars)+propOutliers*dnorm(xOut,mean=meanOutliers,sd=sdOutliers)
outlierDen <- list(x=xOut,y=yOut)

#######################################
# Get the population means
######################################

normalMean <- muNorm
skewMean <- shapeGamma*scaleGamma
superSkewMean <- thetaPareto/(alphaPareto - 1)
outlierMean <- (1-propOutliers)*meanRegulars+propOutliers*meanOutliers

###########################
# Get population standard deviations
##########################

normalSD <- sigmaNorm
skewSD <- sqrt(shapeGamma)*scaleGamma
a <- alphaPareto
superSkewSD <- thetaPareto*sqrt(a/((a-1)^2*(a-2)))
outlierSD <- sqrt((1-propOutliers)*sdRegulars+propOutliers*sdOutliers+(meanOutliers-meanRegulars)^2*propOutliers*(1-propOutliers))

# Define server logic
shinyServer(function(input, output) {
  
  popDen <- reactive({
    switch(input$popDist,
           normal=normalDen,
           skew=skewDen,
           superskew=superSkewDen,
           outliers=outlierDen)
  })
  
  popMean <- reactive({
    switch(input$popDist,
           normal=normalMean,
           skew=skewMean,
           superskew=superSkewMean,
           outliers=outlierMean)
  })
  
  popSD <- reactive({
    switch(input$popDist,
           normal=normalSD,
           skew=skewSD,
           superskew=superSkewSD,
           outliers=outlierSD)
  })
  
  popMax <- reactive({
    switch(input$popDist,
           normal=max(normalDen$x),
           skew=max(skewDen$x),
           superskew=max(superSkewDen$x),
           outliers=max(outlierDen$x)
    )
  })
  
  popMin <- reactive({
    switch(input$popDist,
           normal=min(normalDen$x),
           skew=min(skewDen$x),
           superskew=min(superSkewDen$x),
           outliers=min(outlierDen$y)
    )
  })
  
  yMax <- reactive({
    densities <- popDen()
    max(densities$y)*1.5
  })
  

  simLimit <- 10000 #upper limit on number of sims at once
  
  #Keep track of number of simulations in a given "set-up"
  numberSims <- 0
  xbars <- numeric()
  latestSamp <- numeric()
  latestMean <- numeric()
  
  #we also want the ability to refresh the "set-up
  total <- 0 #total number of sims over all set-ups including current one
  totalPrev <- 0 #total number of sims over all set-ups excluding current one
  
  simsUpdate <- reactive({
    if (input$resample > 0) {


      reps <- min(simLimit,isolate(input$sims))
      total <<- total + reps
      numberSims <<- numberSims + reps
      
      popMean <- isolate(popMean())
      n <- isolate(input$n)
      
      itemNumb <- reps*n
      sampleItems <- isolate(
        switch(input$popDist,
               normal=rnorm(itemNumb,mean=muNorm,sd=sigmaNorm),
               skew=rgamma(itemNumb,shape=shapeGamma,scale=scaleGamma),
               superskew=rpareto(itemNumb,alpha=alphaPareto,theta=thetaPareto),
               outliers=routlier(itemNumb)
        ))
      
      sampleMatrix <- matrix(sampleItems,ncol=n,nrow=reps)
      
      latestSamp <<- sampleMatrix[reps,]
      latestMean <<- mean(latestSamp)
      
      newMeans <- rowSums(sampleMatrix)/n
      
      xbars <<- c(xbars,newMeans)
      
    } #end if resample
      
  })
  
  #this erases the simulation history and puts user back to initial graph
  simsReset <- reactive({
    input$reset
    totalPrev <<- totalPrev + numberSims
    numberSims <<- 0
    xbars <<- numeric()
    latestSamp <<- numeric()
    latestMean <<- numeric()
    return(totalPrev)
  })
  
  #help with conditonal panals
  output$totalPrev <- reactive({
    simsReset()
  })
  
  # needed for the conditional panels to work
  outputOptions(output, 'totalPrev', suspendWhenHidden=FALSE)
  
  output$total <- reactive({
    simsUpdate() #for dependency
    total
  })
  
  # needed for the conditional panels to work
  outputOptions(output, 'total', suspendWhenHidden=FALSE)
  
  
  output$initialGraph <- renderPlot({
    popDen <- popDen()
    popMean <- popMean()
    plot(popDen$x,popDen$y,type="l",lwd=3,col="red",
         main="Density Curve of Population",
         xlab="",
         ylab="density")
    abline(v=popMean,lwd=2)
    
  })
  
  output$xbar <- renderPlot({
    input$resample
    n <- isolate(input$n)
    mu <- isolate(popMean())
    sigma <- isolate(popSD())
    sdMean <- sigma/sqrt(n)
    lower <- mu - 5*sdMean
    upper <- mu + 5*sdMean
    if (numberSims == 1) {
      xbarDen <- density(xbars,n=1024,bw=1)
    }
    if (numberSims >= 2 && n < 5) {
      xbarDen <- density(xbars,n=1024,bw=0.1)
    }
    if (numberSims >= 2 && n >= 5) {
      xbarDen <- density(xbars,n=1024,bw="SJ")
    }
    
    if (numberSims > 0) {
    ymax <- max(xbarDen$y,dnorm(mu,mean=mu,sd=sdMean))
    plot(xbarDen$x,xbarDen$y,type="l",lwd=2,col="blue",
         main="x-bar vs. normal curve",cex.main=2,
         xlab="x-bar", ylim=c(0,ymax),xlim=c(lower,upper),
         ylab="density")
    curve(dnorm(x,mean=mu,sd=sdMean),xlim=c(lower,upper),col="red",lwd=2,add=TRUE)
    if (numberSims <= 50) {
      rug(xbars)
    }
    points(latestMean, 0, col = "blue", pch = 20,cex=2)
    } #end check that there are samples
    
  })
  
  output$graphSample <- renderPlot({
    input$resample #for the dependency
    if(numberSims > 0) {
    popDen <- isolate(popDen())
    popMean <- isolate(popMean())
    n <- isolate(input$n)
    
    samp <- latestSamp
    
    
    xmin <- isolate(popMin())
    xmax <- isolate(popMax())
    ymax <- isolate(yMax())
    
    xbar <- latestMean
    
    hist(samp,freq=FALSE,col="lightblue",xlim=c(xmin,xmax),
         ylim=c(0,ymax),
         main="Population Density Curve\nand Histogram of Sample",
         xlab="x",
         cex.main=2,cex.sub=2
    )
    lines(popDen,lwd=2,col="red")
    abline(v=popMean,lwd=2)
    points(xbar, 0, col = "blue", pch = 20,cex=2)
    
    } # end checking that we have a latest sample
    
  })

  
})
  