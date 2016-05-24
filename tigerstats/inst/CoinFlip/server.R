library(shiny)

source("binomGraphs.R")

# Define server logic for CoinFlip
shinyServer(function(input, output) {
  
  simLimit <- 10000 #upper limit on number of sims at once

  #Keep track of number of simulations in a given "set-up"
  numberSims <- 0
  successSims <- numeric()
  latestSim <- 0
  fullSim <-character()
  
  #we also want the ability to refresh the "set-up
  total <- 0 #total number of sims over all set-ups including current one
  totalPrev <- 0 #total number of sims over all set-ups excluding current one
  
  simsUpdate <- reactive({
    if (input$resample > 0) {
    p <- isolate(input$p)
    n <- isolate(input$n)
    reps <- min(10000,isolate(input$sims))
    newSims <- rbinom(reps,size=n,prob=p)
    successSims <<- c(successSims,newSims)
    latestSim <<- newSims[reps]
    numberSims <<- numberSims + reps
    total <<- total+reps
    
    #now build fake list of outcomes for each trial, on last sim
    full <- c(rep(isolate(input$success),latestSim),
              rep(isolate(input$failure),n-latestSim))
    full <- sample(full,size=n,replace=FALSE)
    fullSim <<- full
    list(numberSims,successSims,latestSim)
    }
  })
  
  #this erases the simulation history and puts user back to initial graph
  simsReset <- reactive({
    input$reset
    totalPrev <<- totalPrev + numberSims
    numberSims <<- 0
    successSims <<- numeric()
    latestSim <<- 0
    return(totalPrev)
  })
  
  
  # for debugging
#   output$successSims <- renderText({
#     simsUpdate()[[2]]
#   })

#print outcome of each trial to screen
output$fullSim <- reactive({
  simsUpdate()
  fullSim
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
  
  latestSim <- reactive({
    simsUpdate()[[3]][1]
  })
  
  valueNames <- reactive({
    c(input$success,input$failure)  
  })
  
  # is our observed value above or below what was expected?
  side <- reactive({
    n <- input$n
    p <- input$p
    x <- input$xObs
    if (x < n*p) return("below") else return("above")
  })
  
  
  output$barGraphInitial <- renderPlot({
    x <- input$xObs
    n <- input$n
    p <- input$p
    
    validate(
      need(input$success,"I'll wait until you enter a name for Success.")
    )
    validate(
      need(input$failure,"I'll wait until you enter a name for Failure.")
    )
    
    validate(
      need(x <= n,"The number of successes should not exceed the number of trials.")
      )
    
    vals <- valueNames()
    observed <- c(x,n-x)
    expected <- c(n*p,n*(1-p))
    tab <- rbind(observed,expected)
    rownames(tab) <-c("Observed","Expected")
    colnames(tab) <- vals
    barplot(tab,beside=T,col=c("#ee7700","grey"),
            main="Bargraph of Observed and Expected Counts",xlab="",ylab="Counts",
            legend.text=TRUE,
            sub=paste("The bar graph above compares the number of successes and failures that were actually\nobserved with the numbers that you would expected\nif the chance of success on each trial were",
                      100*input$p," percent."))
  })
  
  output$remarksInitial <- renderText({
    paste("The bar graph above compares the number of successes and failures that were actually",
          "observed with the numbers that you would expected if the chance of success on each trial were",
          100*input$p," percent.")
  })
  
  output$remarksLatest1 <- renderText({
    simsUpdate() # for the dependency
    paste("Number of successes in this simulation =  ",latestSim,".",sep="")
  })
  
  output$remarksLatest2 <- renderText({
    simsUpdate() # for the dependency
    paste("Number of successes in the most recent simulation =  ",latestSim,".",sep="")
  })
  
  
  output$barGraphLatest <- renderPlot({
    
    input$resample #gets the dependency
    
    x <- isolate(input$xObs)
    n <- isolate(input$n)
    p <- isolate(input$p)
    vals <- isolate(valueNames())
    
    if (numberSims > 0) {
    
    observed <- c(x,n-x)
    expected <- c(n*p,n*(1-p))
    simulated <- c(latestSim,n-latestSim)
    
    tab <- rbind(observed, expected,simulated)
    rownames(tab) <-c("Observed","Expected","Simulated")
    colnames(tab) <- isolate(valueNames())
    barplot(tab,beside=T,col=c("#ee7700","grey","#3333ff"),
            main="Bargraph of Observed, Expected, and Latest Simulation",xlab="",
            ylab="Counts",
            legend.text=TRUE)
  } # end if
  
    
  })
  
  output$histogram <-
    renderPlot({
      input$resample #for the dependency
    minSim <- min(successSims)
    maxSim <- max(successSims)
    myBreaks <- seq(minSim-0.5,maxSim+0.5,by=1)
    histinfo <- hist(successSims,breaks=myBreaks,
       xlab="Number of Successes",
       main="Distribution of Simulated Success-Counts",
       sub=paste0("The red dot shows the results of the last simulation.")
          )
    nvals <- minSim:maxSim
    if (isolate(side()) == "above") {
      Shading <- ifelse(nvals >= isolate(input$xObs),"lightblue",NA)
      } else {
        Shading <- ifelse(nvals <= isolate(input$xObs),"lightblue",NA)
      }
    rect(nvals-0.5,rep(0,times=length(nvals)),nvals+0.5,histinfo$counts,
         col=Shading,border="black")
    points(latestSim,0,col="red",pch=19,cex=2)
    abline(v=isolate(input$xObs))
           
    })
    
output$summary1 <- renderTable({
  simsUpdate() #for the dependency
  observed <- isolate(input$xObs)
  number <- numberSims
  if (side() == "above") {
    extremeCount <- length(successSims[successSims >= observed])
    } else {
    extremeCount <- length(successSims[successSims <= observed])
    }
    extremeCount <- as.integer(extremeCount)
  percent <- paste0(round(100*extremeCount/number,2),"%")
  if (side() == "above") {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulations so far","Number >= Observed","Percentage")
  } else {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulations so far","Number <= observed","Percentage")
  }
  tab
})

output$summary2 <- renderTable({
  simsUpdate() #for dependency
  observed <- isolate(input$xObs)
  number <- numberSims
  if (side() == "above") {
    extremeCount <- length(successSims[successSims >= observed])
  } else {
    extremeCount <- length(successSims[successSims <= observed])
  }
  extremeCount <- as.integer(extremeCount)
  percent <- paste0(round(100*extremeCount/number,2),"%")
  if (side() == "above") {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulations so far","Number >= Observed","Percentage")
  } else {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulations so far","Number <= Observed","Percentage")
  }
  tab
})
  
output$probHistogram <- renderPlot({
  if (side() == "above") {
  binomGraphs(bound=input$xObs-1,region="above",size=input$n,prob=input$p,
              xlab="Number of Successes")
  
  } else {
    binomGraphs(bound=input$xObs,region="below",size=input$n,prob=input$p,
                xlab="Number of Successes")
  }
  abline(v=isolate(input$xObs))
})

output$remarkProb <- renderText({
  if (side() == "above") {
  paste0("The histogram above is what you would get if you could simulate many, many times.",
        " The shaded area gives the exact probability of getting ",
        input$xObs," or more successes in ",input$n," trials, if the chance of success on each trial is ",
        100*input$p," percent.")
  } else {
    paste0("The histogram above is what you would get if you could simulate many, many times.",
           " The area of the shaded recatangles gives the exact probability of getting ",
           input$xObs," or fewer successes in ",input$n," trials, if the chance of success on each trial is ",
           100*input$p," percent.  If you can't see a shaded area, that's because it was too small to show.")
  }
})
  
})
  
