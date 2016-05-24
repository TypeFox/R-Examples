library(shiny)

source("normGraph.R")
source("UnderShade.R")

# Define server logic for CoinFlip
shinyServer(function(input, output) {
  
  simLimit <- 10000 #upper limit on number of sims at once

  #Keep track of number of simulations in a given "set-up"
  numberSims <- 0
  successSims <- numeric()
  allProps <- numeric()
  latestSim <- 0
  
  # get intial info (needed for random \reassignment)
  
  groupSizes <- reactive({
    as.integer(unlist(strsplit(input$groupSizes,split=",")))
  })
  
  successCounts <- reactive({
    as.integer(unlist(strsplit(input$successCounts,split=",")))
  })
  
  groupNames <- reactive({
    unlist(strsplit(input$groupNames,split=","))
  })
  
  goodSizes <- reactive({
    sizes <- as.integer(unlist(strsplit(input$groupSizes,split=",")))
    goodSizes <- TRUE
    if (any(is.na(sizes))) goodSizes <- FALSE
    if (length(sizes) != 2) goodSizes <- FALSE
    if (any(sizes < 1)) goodSizes <- FALSE
    goodSizes
  })
  
  goodCounts <- reactive({
    counts <- as.integer(unlist(strsplit(input$successCounts,split=",")))
    goodCounts <- TRUE
    if (any(is.na(counts))) goodCounts <- FALSE
    if (length(counts) != 2) goodCounts <- FALSE
    if (any(counts < 0)) goodCounts <- FALSE
    goodCounts
  })
  
  goodNames <- reactive({
    names <- unlist(strsplit(input$groupNames,split=","))
    goodNames <- TRUE
    if (any(is.na(names))) goodNames <- FALSE
    if (length(names) != 2) goodNames <- FALSE
    goodNames
  })
  
side <- reactive({
  obsProps <- successCounts()/groupSizes()
  if (obsProps[1] <= obsProps[2]) return("below") else return("above")
})
  
  obsDiff <- reactive({
    groupSizes <- groupSizes()
    n1 <- groupSizes[1]
    n2 <- groupSizes[2]
    successCounts <- successCounts()
    x1 <- successCounts[1]
    x2 <- successCounts[2]
    return(x1/n1-x2/n2)
  })
  
  
  #we also want the ability to refresh the "set-up"
  total <- 0 #total number of sims over all set-ups including current one
  totalPrev <- 0 #total number of sims over all set-ups excluding current one
  
  simsUpdate <- reactive({
    if (input$resample > 0) {
    groupSizes <- isolate(groupSizes())
    n1 <- groupSizes[1]
    n2 <- groupSizes[2]
    successCounts <- isolate(successCounts())
    x1 <- successCounts[1]
    x2 <- successCounts[2]
    
    reps <- min(simLimit,isolate(input$sims))
    
    #think hypergemetric:
    white <- x1+x2
    urn <- n1+n2
    black <- urn - white
    
    newSims <- rhyper(nn=reps,m=white,n=black,k=n1) #number of successes assigned to first group
    successSims <<- c(successSims,newSims)
    allProps <<- successSims/n1-(x1+x2-successSims)/n2
    latestSim <<- newSims[reps]
    numberSims <<- numberSims + reps
    total <<- total+reps
    
    } #end if
    
    })

  
  diffProps <- reactive({
    simsUpdate() # for the dependency
    groupSizes <- isolate(groupSizes())
    n1 <- groupSizes[1]
    n2 <- groupSizes[2]
    successCounts <- isolate(successCounts())
    x1 <- successCounts[1]
    x2 <- successCounts[2]
    return(latestSim/n1-(x1+x2-latestSim)/n2)
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
  
  
  output$barGraphInitial <- renderPlot({
    
    validate(
      need(goodNames(),"")
    )
    
    validate(
      need(goodSizes(),"")
    )
    
    validate(
      need(goodCounts(),"")
    )
    
    validate(
      need(all(successCounts()<=groupSizes()),"")
    )
    
    groupSizes <- groupSizes()
    n1 <- groupSizes[1]
    n2 <- groupSizes[2]
    n <- n1+n2
    successCounts <- successCounts()
    x1 <- successCounts[1]
    x2 <- successCounts[2]
    
    p <- (x1+x2)/n
    
    groupNames <- groupNames()
    group1 <- groupNames[1]
    group2 <- groupNames[2]
    
    success <- input$success
    
    observed <- c(x1,x2)
    expected <- c(p*n1,p*n2)
    tab <- rbind(observed,expected)
    rownames(tab) <-c("Observed","Expected")
    colnames(tab) <- c(paste(group1,success," count"),
                       paste(group2,success," count"))
    barplot(tab,beside=T,col=c("#ee7700","grey"),
            main="Bargraph of Observed and Expected Counts",xlab="",ylab="Counts",
            legend.text=TRUE)
  })
  
  output$remarksInitial <- renderText({
    
    validate(
      need(goodNames(),"")
    )
    
    validate(
      need(goodSizes(),"")
    )
    
    validate(
      need(goodCounts(),"")
    )
    
    validate(
      need(all(successCounts()<=groupSizes()),"")
    )
    
    paste("The bar graph above compares the number of successes actually observed in each group",
          "with the numbers that you would expect if group makes no difference for the response.")
  })

output$remarksInitialMore <- renderText({
  
  validate(
    need(goodNames(),"")
  )
  
  validate(
    need(goodSizes(),"")
  )
  
  validate(
    need(goodCounts(),"")
  )
  
  validate(
    need(all(successCounts()<=groupSizes()),"")
    )
  
    diff <- obsDiff()
    
    groupNames <- isolate(groupNames())
    group1 <- groupNames[1]
    group2 <- groupNames[2]
    
    paste0("In this experiment, the difference in sample proportions (",
           group1," - ",group2,") = ",round(100*diff,2),
           " percent.")
  
})

output$initialTwoWay <- renderTable({
  
  
  validate(
    need(goodNames(),"Enter exactly two group names, separated by a comma.")
  )
  
  validate(
    need(goodSizes(),"Enter the group sizes as two positive whole numbers, separated by a comma.")
  )
  
  validate(
    need(goodCounts(),"Enter the success counts as two non-negative whole numbers, separated by a comma.")
  )
  
  validate(
    need(all(successCounts()<=groupSizes()),
         paste0("One of the groups has more successes than it has members! ",
                "Did you enter the number in the correct order?")
    )
  )
  
  groupSizes <- groupSizes()
  n1 <- groupSizes[1]
  n2 <- groupSizes[2]
  n <- n1+n2
  successCounts <- successCounts()
  x1 <- successCounts[1]
  x2 <- successCounts[2]
  groupNames <- groupNames()
  group1 <- groupNames[1]
  group2 <- groupNames[2]
  
  p1hat <- round(x1/n1,4)
  p2hat <- round(x2/n2,4)
  
  tab <- rbind(c(as.character(x1),as.character(n1-x1),as.character(p1hat)),
               c(as.character(x2),as.character(n2-x2),as.character(p2hat)))
  rownames(tab) <- groupNames
  colnames(tab) <- c(input$success,input$failure,"Proportion")
  tab
})
  
  output$remarksLatest1 <- renderText({

    diff <- diffProps()
    
    groupNames <- isolate(groupNames())
    group1 <- groupNames[1]
    group2 <- groupNames[2]
    
    paste0("For the latest simulation, the difference in sample proportions (",
      group1," - ",group2,") = ",round(100*diff,2),
      " percent.")
    

  })

output$latestTwoWayBar <- renderTable({
  simsUpdate() #for the dependency
  groupSizes <- isolate(groupSizes())
  n1 <- groupSizes[1]
  n2 <- groupSizes[2]
  successCounts <- isolate(successCounts())
  x1 <- successCounts[1]
  x2 <- successCounts[2]
  groupNames <- isolate(groupNames())
  group1 <- groupNames[1]
  group2 <- groupNames[2]
  
  p1hat <- round(latestSim/n1,4)
  p2hat <- round((x1+x2-latestSim)/n2,4)
  
  tab <- rbind(c(as.character(latestSim),as.character(n1-latestSim),as.character(p1hat)),
               c(as.character(x1+x2-latestSim),as.character(n2-x1-x2+latestSim),as.character(p2hat)))
  rownames(tab) <- groupNames
  colnames(tab) <- c(input$success,input$failure,"Proportion")
  tab
})
  
  output$remarksLatest2 <- renderText({
    

    diff <- diffProps()
    
    groupNames <- isolate(groupNames())
    group1 <- groupNames[1]
    group2 <- groupNames[2]
    
    paste0("For the latest simulation, the difference in sample proportions (",
           group1," - ",group2,") = ",round(100*diff,2),
                                             " percent.")
    

    
  })

output$latestTwoWayDen <- renderTable({
  simsUpdate() #for the dependency
  groupSizes <- isolate(groupSizes())
  n1 <- groupSizes[1]
  n2 <- groupSizes[2]
  successCounts <- isolate(successCounts())
  x1 <- successCounts[1]
  x2 <- successCounts[2]
  groupNames <- isolate(groupNames())
  group1 <- groupNames[1]
  group2 <- groupNames[2]
  
  p1hat <- round(latestSim/n1,4)
  p2hat <- round((x1+x2-latestSim)/n2,4)
  
  tab <- rbind(c(as.character(latestSim),as.character(n1-latestSim),as.character(p1hat)),
               c(as.character(x1+x2-latestSim),as.character(n2-x1-x2+latestSim),as.character(p2hat)))
  rownames(tab) <- groupNames
  colnames(tab) <- c(input$success,input$failure,"Proportion")
  tab
})
  
  
  output$barGraphLatest <- renderPlot({
    
    input$resample #gets the dependency
    
    groupSizes <- isolate(groupSizes())
    n1 <- groupSizes[1]
    n2 <- groupSizes[2]
    n <- n1+n2
    successCounts <- isolate(successCounts())
    x1 <- successCounts[1]
    x2 <- successCounts[2]
    
    groupNames <- isolate(groupNames())
    group1 <- groupNames[1]
    group2 <- groupNames[2]
    
    success <- isolate(input$success)
    
    if (numberSims > 0) {
    
    p <- (x1+x2)/n
      
    observed <- c(x1,x2)
    expected <- c(p*n1,p*n2)
    simulated <- c(latestSim,x1+x2-latestSim)
    tab <- rbind(observed,expected,simulated)  
      
    rownames(tab) <-c("Observed","Expected","Simulated")
    colnames(tab) <- c(paste0(success," count (",group1,")"),
                       paste0(success," count (",group2,")"))
    barplot(tab,beside=T,col=c("#ee7700","grey","#3333ff"),
            main="Bargraph of Observed, Expected, and Latest Simulation",xlab="",
            ylab="Counts",
            legend.text=TRUE)
  } # end if
  
    
  })
  
  output$density <-
  
# if we want to go back to using a historgram:
#     renderPlot({
#     diff <- diffProps() #also gives us the dependency
#     
#     hist(allProps,
#        xlab="Difference",
#        main="Distribution of Differences in Group Proportions"
#           )
#     points(diff,0,col="red",pch=19,cex=2)
#     
#     abline(v=isolate(obsDiff()))
    
    renderPlot({
      input$resample
      if (length(allProps)==1) band <- 1 else band <- "nrd0"
      dprops <- density(allProps,n=500,bw=band)
      plot(dprops$x,dprops$y,type="l",col="blue",
           xlab="Difference",ylab="Estimated Density",
           main="Distribution of Resampled Differences")
      latest <- allProps[length(allProps)]
      points(latest,0,col="red",pch=19)
      abline(v=isolate(obsDiff()))
           
    })
    
output$summary1 <- renderTable({
  simsUpdate() #for the dependency

  observed <- isolate(obsDiff())
  number <- numberSims
  if (side() == "above") {
    extremeCount <- length(allProps[allProps >= observed])
    } else {
    extremeCount <- length(allProps[allProps <= observed])
    }
    extremeCount <- as.integer(extremeCount)
  percent <- paste0(round(100*extremeCount/number,2),"%")
  if (side() == "above") {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulated Differences","Number >= Observed Diff","Percentage")
  } else {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulations","Number <= Observed Diff","Percentage")
  }
  tab
  

})

output$summary2 <- renderTable({
  simsUpdate() #for the dependency
  

  observed <- isolate(obsDiff())
  
  number <- numberSims
  if (side() == "above") {
    extremeCount <- length(allProps[allProps >= observed])
  } else {
    extremeCount <- length(allProps[allProps <= observed])
  }
  extremeCount <- as.integer(extremeCount)
  percent <- paste0(round(100*extremeCount/number,2),"%")
  if (side() == "above") {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulated Differences","Number >= Observed Diff","Percentage")
  } else {
    tab <- data.frame(as.integer(number),extremeCount,percent)
    names(tab) <- c("Simulations","Number <= Observed Diff","Percentage")
  }
  tab
  

})
 
output$normalCurve <- renderPlot({

  groupSizes <- isolate(groupSizes())
  n1 <- groupSizes[1]
  n2 <- groupSizes[2]
  successCounts <- isolate(successCounts())
  x1 <- successCounts[1]
  x2 <- successCounts[2]
  
  observed <- isolate(obsDiff())
  
  #This needs to be improved (or should it be??)
  p <- (x1+x2)/(n1+n2)
  sdDiff <- sqrt(p*(1-p)*(1/n1+1/n2))
  
  
  if (side() == "above") {
  normGraph(bound=observed,region="above",mean=0,sd=sdDiff)
  
  } else {
    normGraph(bound=observed,region="below",mean=0,sd=sdDiff)
  }
  

})

output$remarksProb <- renderText({
  
  observed <- obsDiff()
  if (side() == "above") {
  paste0("The curve above approximates the histogram you would get if you could simulate many, many times.",
        " The shaded area gives an approximate probability of getting an observed difference in proportions",
       "of ",round(observed,3)," or more, if group makes no difference in the response.")
  } else {
    paste0("The curve above approximates the histogram you would get if you could simulate many, many times.",
           " The shaded area gives an approximate probability of getting an observed difference in proportions",
           "of ",round(observed,3)," or less, if group makes no difference in the response.")
  }
})
  
})
  
