library(shiny)
library(shape)

N <- 52 #sqrt of pop size

xvals <- rep(1:N,each=N)
yvals <- rep(1:N,times=N)
xreg <- floor((xvals-1)/N*4)
yreg <- floor((yvals-1)/N*4)

stratum <- character(N^2)
for (i in 1:N^2) {
  stratum[i] <- paste0(xreg[i],",",yreg[i])
}

spots <- unique(stratum)

status <- rep(FALSE,N^2)
pointpop <- data.frame(x=xvals,y=yvals,stratum=stratum,status=status)




# Define server logic for SamplingMethods
shinyServer(function(input, output) {

  
  
  #we also want the ability to refresh the "set-up"
  total <- 0 #total number of samples over all set-ups including current one
  totalPrev <- 0 #total number of sanples over all set-ups excluding current one
  currentCount <- 0
  

  sampUpdate <- reactive({
    if (input$sample > 0) {

    currentCount <<- currentCount +1
    total <<- total + 1
    }
    
  })
  
  
  #this erases the history and puts user back to initial graph
  appReset <- reactive({
    input$reset
    totalPrev <<- totalPrev + currentCount
    currentCount <<- 0
    return(totalPrev)
  })


  
  #help with conditonal panals
  output$totalPrev <- reactive({
    appReset()
  })
  
  # needed for the conditional panels to work
  outputOptions(output, 'totalPrev', suspendWhenHidden=FALSE)
  
  output$total <- reactive({
    sampUpdate() #for dependency
    total
  })
  
  # needed for the conditional panels to work
  outputOptions(output, 'total', suspendWhenHidden=FALSE)
  
  output$population <- renderPlot({
    input$n #to get us going
    
    plot(0,0,axes=F,xlim=c(-0.7,N+0.7),ylim=c(-0.7,N+0.7),col="transparent",
         xlab="",ylab="",main="Say Hello to the Population!")
    rect(0,0,N+0.7,N+0.7,lwd=2)
    
    if (!input$hide) {
      for (i in 1:3) {
        lines(x=c((N+1)*i/4,(N+1)*i/4),y=c(0,N+0.7))
        lines(x=c(0,N+0.7),y=c((N+1)*i/4,(N+1)*i/4))
      }
    }
    
    points(pointpop$x,pointpop$y,pch=19,cex=0.1)
    
    
    
  })
  
  
  output$sampleplot <- renderPlot({

    sampUpdate()
    print(total)
    method <- input$method
    n <- as.numeric(input$n)
    
    meth <- switch(method,
                   "srs" = "Simple Random",
                   "strat" = "Stratified",
                   "cluster" = "Cluster"
      )
    
    plot(0,0,axes=F,xlim=c(-0.7,N+0.7),ylim=c(-0.7,N+0.7),col="transparent",
         xlab="",ylab="",main=paste0("A ",meth," Sample"))
    
    rect(0,0,N+0.7,N+0.7,lwd=2)
    
    if (!input$hide) {
      for (i in 1:3) {
        lines(x=c((N+1)*i/4,(N+1)*i/4),y=c(0,N+0.7))
        lines(x=c(0,N+0.7),y=c((N+1)*i/4,(N+1)*i/4))
      }
    }
      
    if (method=="srs") {
    
      centers <- pointpop[sample(1:(N^2),size=n,replace=F),]
      for (i in 1:n) {
        midpoint <- c(centers$x[i],centers$y[i])
        print(filledcircle(r1=0.25,r2=0.0,mid=midpoint,col="blue"))
      }
    }
    
    
    if (method=="strat") {
      
      m <- N^2/16
      chosen <- c(rep(TRUE,n/16),rep(FALSE,m-n/16))
      
      
      for (spot in spots) {
        spotPicks <- sample(chosen,size=m,replace=FALSE)
        pointpop$status[pointpop$stratum==spot] <- spotPicks
      }
      
      centers <- pointpop[pointpop$status,]
      for (i in 1:n) {
        midpoint <- c(centers$x[i],centers$y[i])
        print(filledcircle(r1=0.25,r2=0.0,mid=midpoint,col="blue"))
      }
    }
    
    
    if (method=="cluster") {
      
      #pick your four clusters
      selectSpots <- sample(spots,size=4,replace=FALSE)
      
      tot <- N^2/16 #total in a cluster
      pickNum <- n/4 #number to pick from each cluster
      chosen <- c(rep(TRUE,pickNum),rep(FALSE,tot - pickNum))
      
      for (spot in selectSpots) {
        spotPicks <- sample(chosen,size=length(chosen),replace=FALSE)
        pointpop$status[pointpop$stratum==spot] <- spotPicks
      }
      
      centers <- pointpop[pointpop$status,]
      for (i in 1:n) {
        midpoint <- c(centers$x[i],centers$y[i])
        print(filledcircle(r1=0.25,r2=0.0,mid=midpoint,col="blue"))
      }
    }
    
  })
  
  output$explanation <- renderText({
    
    method <- input$method
    n <- as.numeric(input$n)
    
    if (method=="srs") {
      return(paste0("\n\nEach set of ",n," items in the population has the same chance to be the",
             " sample selected.  Items are selected without regard to their districts."))
    }
    
    if (method=="strat") {
      return(paste0("\n\nFor each district, we sampled ",n/16," item(s) at random."))
    }
    
    if (method=="cluster") {
      return(paste0("\n\nWe picked four districts at random.  From each selected district, we chose ",
             n/4," items at random, making four 'clusters'."))
    }
    
  })
  
})
  
