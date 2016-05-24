### server.R
source("miscFunctions.R")

names <- privateSales$Country

### Read in WHO death data
WHOdeaths <- u5deathsWHO2012
WHOmalaria <- u5malariadeathsWHO2010

# get colors
palette <- RColorBrewer::brewer.pal(9, "Blues")

#print("I Am HERE FIRST")

# drop row with totals
privateSales <- privateSales[-length(privateSales$Country), ]

# rename columns
colnames(privateSales) <- c("Country", "pfPositiveSales", "StdDev")
numCountries <- length(privateSales$Country)

# create cohen data table to be printed on website
displaySales <- privateSales
colnames(displaySales) <- c("Country", "Mean", "Standard_Deviation")
displaySales$Mean <- as.integer(displaySales$Mean)
displaySales$Standard_Deviation <- as.integer(displaySales$Standard_Deviation)

## global random 6 digit number that is generated for each set of input parameters
inputID <- 100000

### Define server logic for slider examples ###
shinyServer(function(input, output) {

  #print("I am here Second")

  generateLHS <- reactive({
    input$goRefresh # make dependent on Refresh Button
    inputID <- runif(1, 0, 999999)
    inputID <- as.integer(inputID)
    #print(paste(inputID, "sim1"))
    #print(inputID)
    sameFlag <- FALSE

    countrySpec <- input$countrySpecific
    manuscript <- input$manuscript
    numReps <- input$N
    minCFR <- input$CFR[1]
    maxCFR <- input$CFR[2]
    if(minCFR == maxCFR) sameFlag <- TRUE
    surveyCountries = c(29,33,32,35,25,39,15,34,37)

    if(countrySpec){
      # make dependent on Update button
      input$goButton

      # use isolate function, so output doesn't refresh until you click update button
      minPrev <- isolate(c(input$Prev1[1],input$Prev2[1],input$Prev3[1],input$Prev4[1],input$Prev5[1],input$Prev6[1],input$Prev7[1],input$Prev8[1],input$Prev9[1],input$Prev10[1],input$Prev11[1],input$Prev12[1],input$Prev13[1],input$Prev14[1],input$Prev15[1],input$Prev16[1],input$Prev17[1],input$Prev18[1],input$Prev19[1],input$Prev20[1],input$Prev21[1],input$Prev22[1],input$Prev23[1],input$Prev24[1],input$Prev25[1],input$Prev26[1],input$Prev27[1],input$Prev28[1],input$Prev29[1],input$Prev30[1],input$Prev31[1],input$Prev32[1],input$Prev33[1],input$Prev34[1],input$Prev35[1],input$Prev36[1],input$Prev37[1],input$Prev38[1],input$Prev39[1]))

      maxPrev <- isolate(c(input$Prev1[2],input$Prev2[2],input$Prev3[2],input$Prev4[2],input$Prev5[2],input$Prev6[2],input$Prev7[2],input$Prev8[2],input$Prev9[2],input$Prev10[2],input$Prev11[2],input$Prev12[2],input$Prev13[2],input$Prev14[2],input$Prev15[2],input$Prev16[2],input$Prev17[2],input$Prev18[2],input$Prev19[2],input$Prev20[2],input$Prev21[2],input$Prev22[2],input$Prev23[2],input$Prev24[2],input$Prev25[2],input$Prev26[2],input$Prev27[2],input$Prev28[2],input$Prev29[2],input$Prev30[2],input$Prev31[2],input$Prev32[2],input$Prev33[2],input$Prev34[2],input$Prev35[2],input$Prev36[2],input$Prev37[2],input$Prev38[2],input$Prev39[2]))

      for(i in 1:39){
        if(minPrev[i] == maxPrev[i]) sameFlag <- TRUE
      }
    }
    else{
      minPrev <- input$Prev[1]
      maxPrev <- input$Prev[2]
      if(minPrev == maxPrev) sameFlag <- TRUE
    }


    ## Get sample vector for each variable
    deathRate = getUniformLHS(numReps, minCFR, maxCFR) # 1 in 1000 - 5 in 1000
    sales = matrix(NA, nrow = numReps, ncol = numCountries)
    fakePercent = matrix(NA, nrow = numReps, ncol = numCountries)

    if(countrySpec){
      for (i in 1:numCountries) {
        sales[, i] = getNormalLHS(numReps, privateSales$pfPositiveSales[i], privateSales$StdDev[i])
        fakePercent[, i] = getUniformLHS(numReps, minPrev[i], maxPrev[i])
      }
    } else if(manuscript){
      j = 1
      for (i in 1:numCountries) {
        sales[, i] = getNormalLHS(numReps, privateSales$pfPositiveSales[i], privateSales$StdDev[i])

        if (i %in% surveyCountries){
          fakePercent[, i] = getNormalLHS(numReps, antimalarialQualitySurvey$Mean[j], antimalarialQualitySurvey$SE[j], belowOneOnly = TRUE)
          j = j + 1
        } else {
          fakePercent[, i] = getUniformLHS(numReps, 0, 0.40)
        }
      }
    } else {
      j = 1
      for (i in 1:numCountries) {
        sales[, i] = getNormalLHS(numReps, privateSales$pfPositiveSales[i], privateSales$StdDev[i])
        fakePercent[, i] = getUniformLHS(numReps, minPrev, maxPrev)
      }
    }

    calculations <- matrix(NA, nrow = numReps, ncol = (numCountries + 1))

    ## Run simulations
    for (i in 1:numReps) {
      calculations[i, ] = model(sales[i, ], fakePercent[i, ], deathRate[i])

      for (j in 1:40) {
        if (calculations[i, j] < 0.5) {
          calculations[i, j] = 0
        }
      }
    }

    #calculations
    calculations <- data.frame(calculations)
    colnames(calculations) <- names
    calculations <- list(calculations, deathRate, sales, fakePercent,inputID, sameFlag)
    calculations
  }
  )

  generateSummary <- reactive({
    numReps <- input$N

    results <- generateLHS()
    results <- results[1]; results <- data.frame(results)

    #results <- data.frame(results)
    totalDeaths <- results[,40]
    results[, 41] = 1:numReps
    colnames(results) = c(paste(privateSales$Country), "All Countries", "SimulationNumber")

    # reformat data from wide to long format
    meltedOutput <- reshape2::melt(results, id.vars = c("SimulationNumber"))

    # use summary function to generate summary
    meltedSummary <- summarySE(meltedOutput, measurevar = "value", groupvars = c("variable"))
    meltedSummary <- meltedSummary[order(meltedSummary$value, decreasing = TRUE),]

    #print(colnames(meltedSummary))
    #print(meltedSummary)

    meltedSummary <- meltedSummary[,c(1,3,4,5,6,7,8,9)]

    meltedSummary <- meltedSummary[,c(1,5,6,8,2,7,4,3)]

    meltedSummary[,2] <- as.integer(meltedSummary[,2])
    meltedSummary[,3] <- as.integer(meltedSummary[,3])
    meltedSummary[,4] <- as.integer(meltedSummary[,4])
    meltedSummary[,5] <- as.integer(meltedSummary[,5])
    meltedSummary[,6] <- as.integer(meltedSummary[,6])
    meltedSummary[,7] <- as.integer(meltedSummary[,7])
    meltedSummary[,8] <- as.integer(meltedSummary[,8])

    colnames(meltedSummary) <- c("Country", "Min", "First Quartile", "Median","Mean", "Third Quartile", "Max", "Std Dev")
    rownames(meltedSummary) <- NULL

    meltedSummary$Country <- paste(meltedSummary$Country)
    #print(class(meltedSummary$Country))

    for (i in 1:length(meltedSummary$Country)){
      if (meltedSummary$Country[i] == "Cote d'Ivoire"){
        meltedSummary$Country[i] <- "Coté d'Ivoire"
      }
    }

    #write.csv(meltedSummary$Country,"testNames.csv")
    #test()

    meltedSummary

  }
  )

  # generate an HTML table view of the summary data
  output$summaryTable <- renderTable({
    data.frame(generateSummary())
  }
  )

  # function to allow for .csv download of summary stats
  output$downloadSummary <- downloadHandler(
    filename = function() {
      results <- generateLHS()
      inputID <- results[5]
      #print(inputID)
      paste(inputID,'_SummaryStats.csv', sep='')
    }
    ,content = function(file){
      out <- generateSummary()
      for (i in 1:length(out[,1])){
        if(out[i,1] == "Coté d'Ivoire") out[i,1] <- "Cote d'Ivoire"
      }
      write.csv(out, file)
    }
  )

  plotRaw <- reactive({
    numReps <- input$N
    results <- generateLHS()
    results <- results[1]; results <- data.frame(results)

    totalDeaths <- results[,40]
    results[, 41] = 1:numReps
    colnames(results) = c(paste(privateSales$Country), "All Countries", "SimulationNumber")

    # reformat data from wide to long format
    meltedOutput <- reshape2::melt(results, id.vars = c("SimulationNumber"))

    # use summary function to generate summary
    meltedSummary <- summarySE(meltedOutput, measurevar = "value", groupvars = c("variable"))

    medianGraphs <- meltedSummary
    medianGraphs2 <- medianGraphs[1:39,]

    # turn names into text
    medianGraphs2$variable <- paste(medianGraphs2$variable)

    # order meltedSummary 2 by country name
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$variable),]

    #print(medianGraphs2$variable)

    # shorten long country names so that they appear on the graph
    for (i in 1:length(colnames(medianGraphs2))){
      if (medianGraphs2$variable[i] == "Central African Republic") {
        medianGraphs2$variable[i] <- "CAR"
      }
      if (medianGraphs2$variable[i] == "Democratic Republic of the Congo") {
        medianGraphs2$variable[i] <- "DRC"
      }
      if (medianGraphs2$variable[i] == "Cote d'Ivoire"){
        medianGraphs2$variable[i] <- "Coté d'Ivoire"
      }
    }

    #print(medianGraphs2$variable)

    # sort by estimated deaths
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$median, decreasing = TRUE),]

    if(medianGraphs2$variable[39] != 'Nigeria') {}

    # get max without Nigeria
    #ymax = max(medianGraphs2$q3[1:38])

    # get max
    ymax = max(medianGraphs2$q3[1:39])

    #print(ymax)
    par(oma = c(3,1,0,0), srt = -40, xpd = NA)
    plot <-barplot(medianGraphs2$median[1:39],xaxs = "i", xlab = NULL, ylab = NULL,xlim = xLimsBarPlot(39,padding = 0.012), ylim = c(0,ymax*1.05), col = NA, border = NA, las = 1, axes = FALSE)
    mtext(side = 2, text = "Estimated Under-Five Deaths Per Year", outer = TRUE)
    text(x = as.vector(plot)-.3, y = -(.04*ymax), medianGraphs2$variable[1:39], adj = c(0,0), cex = .95)
    axis(2, tck=1, ,col.ticks="light gray", labels = FALSE, lty = 2, lwd = 0.75)
    barplot(add = TRUE, medianGraphs2$median[1:39],xaxs = "i",xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = palette[7], las = 1, main = "Median Death Estimates\n(error bars as Interquartile Range)")
    suppressWarnings(arrows(plot, medianGraphs2$q1[1:39], plot, medianGraphs2$q3[1:39], angle = 90, code = 3, length = 0.05))
  }
  )

  plotRaw2 = function(){
    numReps <- input$N
    results <- generateLHS()
    results <- results[1]; results <- data.frame(results)

    totalDeaths <- results[,40]
    results[, 41] = 1:numReps
    colnames(results) = c(paste(privateSales$Country), "All Countries", "SimulationNumber")

    # reformat data from wide to long format
    meltedOutput <-reshape2::melt(results, id.vars = c("SimulationNumber"))

    # use summary function to generate summary
    meltedSummary <- summarySE(meltedOutput, measurevar = "value", groupvars = c("variable"))

    medianGraphs <- meltedSummary
    medianGraphs2 <- medianGraphs[1:39,]

    # turn names into text
    medianGraphs2$variable <- paste(medianGraphs2$variable)

    # order meltedSummary 2 by country name
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$variable),]

    #print(medianGraphs2$variable)

    # shorten long country names so that they appear on the graph
    for (i in 1:length(colnames(medianGraphs2))){
      if (medianGraphs2$variable[i] == "Central African Republic") {
        medianGraphs2$variable[i] <- "CAR"
      }
      if (medianGraphs2$variable[i] == "Democratic Republic of the Congo") {
        medianGraphs2$variable[i] <- "DRC"
      }
      if (medianGraphs2$variable[i] == "Cote d'Ivoire"){
        medianGraphs2$variable[i] <- "Coté d'Ivoire"
      }
    }

    #print(medianGraphs2$variable)

    # sort by estimated deaths
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$median, decreasing = TRUE),]

    if(medianGraphs2$variable[39] != 'Nigeria') {}

    # get max without Nigeria
    #ymax = max(medianGraphs2$q3[1:38])

    # get max
    ymax = max(medianGraphs2$q3[1:39])

    #print(ymax)
    par(oma = c(3,1,0,0), srt = -40, xpd = NA)
    plot <-barplot(medianGraphs2$median[1:39],xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = NA, border = NA, las = 1, axes = FALSE)
    mtext(side = 2, text = "Estimated Under-Five Deaths Per Year", outer = TRUE)
    text(x = as.vector(plot)-.3, y = -(.02*ymax), medianGraphs2$variable[1:39], adj = c(0,0), cex = .85)
    axis(2, tck=1, ,col.ticks="light gray", labels = FALSE, lty = 2, lwd = 0.75)
    barplot(add = TRUE, medianGraphs2$median[1:39],xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = palette[7], las = 1, main = "Median Death Estimates\n(error bars as Interquartile Range)")
    suppressWarnings(arrows(plot, medianGraphs2$q1[1:39], plot, medianGraphs2$q3[1:39], angle = 90, code = 3, length = 0.05))

  }

  output$raw <-  renderPlot({
    plotRaw()
  }
  )

  # handler for downloading median Estimates as pdf
  output$downloadMedianEstimates <- downloadHandler(

    filename = function() {
      results <- generateLHS()
      inputID <- results[5]
      paste(inputID,'_MedianDeathEstimates.pdf', sep='')
    }

    ,content = function(file)
    {
      pdf(file, width = 11, height = 8.5)
      plotRaw2()
      dev.off()
    }

    ,contentType = 'application/pdf'
  )

  plotMalariaProp <- reactive({
    numReps <- input$N
    results <- generateLHS()
    results <- results[1]; results <- data.frame(results)

    totalDeaths <- results[,40]
    results[, 41] = 1:numReps
    colnames(results) = c(paste(privateSales$Country), "All Countries", "SimulationNumber")

    # reformat data from wide to long format
    meltedOutput <-reshape2::melt(results, id.vars = c("SimulationNumber"))

    # use summary function to generate summary
    meltedSummary <- summarySE(meltedOutput, measurevar = "value", groupvars = c("variable"))

    medianGraphs <- meltedSummary
    medianGraphs2 <- medianGraphs[1:39,]

    # turn names into text
    medianGraphs2$variable <- paste(medianGraphs2$variable)

    # order meltedSummary 2 by country name and load in <5 Deaths
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$variable),]
    medianGraphs2["whoDeaths"] <- WHOdeaths[,3]
    medianGraphs2["whoDeathRate"] <- WHOdeaths[,4]

    # load in <5 malaria deaths
    WHOmalaria <- WHOmalaria[WHOmalaria$Country %in% medianGraphs2$variable,]
    WHOmalaria <- WHOmalaria[order(WHOmalaria$Country),]
    medianGraphs2["malariaDeaths"] <- WHOmalaria$Numeric

    # calculation prop of deaths
    medianGraphs2$whoDeaths <- medianGraphs2$whoDeaths * 1000
    medianGraphs2["malariaProp"] <- medianGraphs2$median / medianGraphs2$malariaDeaths
    medianGraphs2["deathProp"] <- medianGraphs2$median / medianGraphs2$whoDeaths

    #print(medianGraphs2$variable)

    # shorten long country names so that they appear on the graph
    for (i in 1:length(colnames(medianGraphs2))){
      if (medianGraphs2$variable[i] == "Central African Republic") {
        medianGraphs2$variable[i] <- "CAR"
      }
      if (medianGraphs2$variable[i] == "Democratic Republic of the Congo") {
        medianGraphs2$variable[i] <- "DRC"
      }
      if (medianGraphs2$variable[i] == "Cote d'Ivoire"){
        medianGraphs2$variable[i] <- "Coté d'Ivoire"
      }
    }

    #print(medianGraphs2$variable)

    # sort by malaria prop
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$malariaProp, decreasing = TRUE),]

    ymax = max(medianGraphs2$q3/medianGraphs2$malariaDeaths)

    #print(ymax)

    par(oma = c(3,1,0,0), srt = -40, xpd = NA)
    plot <-barplot(medianGraphs2$malariaProp,xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = NA, border = NA, las = 1, axes = FALSE)
    mtext(side = 2, text = "Proportion of Under-Five Malaria Deaths", outer = TRUE)
    text(x = as.vector(plot)-.3, y = -(.04*ymax), medianGraphs2$variable, adj = c(0,0), cex = .95)
    axis(2, tck=1, ,col.ticks="light gray", labels = FALSE, lty = 2, lwd = 0.75)
    barplot(add = TRUE, medianGraphs2$malariaProp,xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = palette[7], las = 1, main = "Median Death Estimates as a Proportion\nof Total Under-Five Malaria Deaths\n(error bars as Interquartile Range)")
    suppressWarnings(arrows(plot, medianGraphs2$q1/medianGraphs2$malariaDeaths, plot, medianGraphs2$q3/medianGraphs2$malariaDeaths, angle = 90, code = 3, length = 0.05))

  }
  )

  plotMalariaProp2 = function(){
    numReps <- input$N
    results <- generateLHS()
    results <- results[1]; results <- data.frame(results)

    totalDeaths <- results[,40]
    results[, 41] = 1:numReps
    colnames(results) = c(paste(privateSales$Country), "All Countries", "SimulationNumber")

    # reformat data from wide to long format
    meltedOutput <-reshape2::melt(results, id.vars = c("SimulationNumber"))

    # use summary function to generate summary
    meltedSummary <- summarySE(meltedOutput, measurevar = "value", groupvars = c("variable"))

    medianGraphs <- meltedSummary
    medianGraphs2 <- medianGraphs[1:39,]

    # turn names into text
    medianGraphs2$variable <- paste(medianGraphs2$variable)

    # order meltedSummary 2 by country name and load in <5 Deaths
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$variable),]
    medianGraphs2["whoDeaths"] <- WHOdeaths[,3]
    medianGraphs2["whoDeathRate"] <- WHOdeaths[,4]

    # load in <5 malaria deaths
    WHOmalaria <- WHOmalaria[WHOmalaria$Country %in% medianGraphs2$variable,]
    WHOmalaria <- WHOmalaria[order(WHOmalaria$Country),]
    medianGraphs2["malariaDeaths"] <- WHOmalaria$Numeric

    # calculation prop of deaths
    medianGraphs2$whoDeaths <- medianGraphs2$whoDeaths * 1000
    medianGraphs2["malariaProp"] <- medianGraphs2$median / medianGraphs2$malariaDeaths
    medianGraphs2["deathProp"] <- medianGraphs2$median / medianGraphs2$whoDeaths

    #print(medianGraphs2$variable)

    # shorten long country names so that they appear on the graph
    for (i in 1:length(colnames(medianGraphs2))){
      if (medianGraphs2$variable[i] == "Central African Republic") {
        medianGraphs2$variable[i] <- "CAR"
      }
      if (medianGraphs2$variable[i] == "Democratic Republic of the Congo") {
        medianGraphs2$variable[i] <- "DRC"
      }
      if (medianGraphs2$variable[i] == "Cote d'Ivoire"){
        medianGraphs2$variable[i] <- "Coté d'Ivoire"
      }
    }

    #print(medianGraphs2$variable)

    # sort by malaria prop
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$malariaProp, decreasing = TRUE),]

    ymax = max(medianGraphs2$q3/medianGraphs2$malariaDeaths)

    #print(ymax)

    par(oma = c(3,1,0,0), srt = -40, xpd = NA)
    plot <-barplot(medianGraphs2$malariaProp,xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = NA, border = NA, las = 1, axes = FALSE)
    mtext(side = 2, text = "Proportion of Under-Five Malaria Deaths", outer = TRUE)
    text(x = as.vector(plot)-.3, y = -(.02*ymax), medianGraphs2$variable, adj = c(0,0), cex = .85)
    axis(2, tck=1, ,col.ticks="light gray", labels = FALSE, lty = 2, lwd = 0.75)
    barplot(add = TRUE, medianGraphs2$malariaProp,xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = palette[7], las = 1, main = "Median Death Estimates as a Proportion\nof Total Under-Five Malaria Deaths\n(error bars as Interquartile Range)")
    suppressWarnings(arrows(plot, medianGraphs2$q1/medianGraphs2$malariaDeaths, plot, medianGraphs2$q3/medianGraphs2$malariaDeaths, angle = 90, code = 3, length = 0.05))

  }

  output$malariaProp <- renderPlot({
    plotMalariaProp()
  }
  )

  # handler for downloading malariaProp as pdf
  output$downloadMalariaProp <- downloadHandler(

    filename = function() {
      results <- generateLHS()
      inputID <- results[5]
      #print(inputID)
      paste(inputID,'_ProportionOfTotalMalariaDeaths.pdf', sep='')
      #name <- "HistogramTest.png"
      #name
    }

    ,content = function(file)
    {
      pdf(file, width = 11, height = 8.5)
      plotMalariaProp2()
      #plot(1,1)
      #hist(runif(1000, 1, 10))
      dev.off()
    }

    ,contentType = 'application/pdf'
  )

  plotDeathProp <-reactive({
    numReps <- input$N
    results <- generateLHS()
    results <- results[1]; results <- data.frame(results)

    totalDeaths <- results[,40]
    results[, 41] = 1:numReps
    colnames(results) = c(paste(privateSales$Country), "All Countries", "SimulationNumber")

    # reformat data from wide to long format
    meltedOutput <-reshape2::melt(results, id.vars = c("SimulationNumber"))

    # use summary function to generate summary
    meltedSummary <- summarySE(meltedOutput, measurevar = "value", groupvars = c("variable"))

    medianGraphs <- meltedSummary
    medianGraphs2 <- medianGraphs[1:39,]

    # turn names into text
    medianGraphs2$variable <- paste(medianGraphs2$variable)

    # order meltedSummary 2 by country name and load in <5 Deaths
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$variable),]
    medianGraphs2["whoDeaths"] <- WHOdeaths[,3]
    medianGraphs2["whoDeathRate"] <- WHOdeaths[,4]

    # load in <5 malaria deaths
    WHOmalaria <- WHOmalaria[WHOmalaria$Country %in% medianGraphs2$variable,]
    WHOmalaria <- WHOmalaria[order(WHOmalaria$Country),]
    medianGraphs2["malariaDeaths"] <- WHOmalaria$Numeric

    # calculation prop of deaths
    medianGraphs2$whoDeaths <- medianGraphs2$whoDeaths * 1000
    medianGraphs2["malariaProp"] <- medianGraphs2$median / medianGraphs2$malariaDeaths
    medianGraphs2["deathProp"] <- medianGraphs2$median / medianGraphs2$whoDeaths

    #print(medianGraphs2$variable)

    # shorten long country names so that they appear on the graph
    for (i in 1:length(colnames(medianGraphs2))){
      if (medianGraphs2$variable[i] == "Central African Republic") {
        medianGraphs2$variable[i] <- "CAR"
      }
      if (medianGraphs2$variable[i] == "Democratic Republic of the Congo") {
        medianGraphs2$variable[i] <- "DRC"
      }
      if (medianGraphs2$variable[i] == "Cote d'Ivoire"){
        medianGraphs2$variable[i] <- "Coté d'Ivoire"
      }
    }

    #print(medianGraphs2$variable)

    # sort by death prop
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$deathProp, decreasing = TRUE),]

    ymax = max(medianGraphs2$q3/medianGraphs2$whoDeaths)

    #print(ymax)

    par(oma = c(3,1,0,0), srt = -40, xpd = NA)
    plot <-barplot(medianGraphs2$deathProp,xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = NA, border = NA, las = 1, axes = FALSE)
    mtext(side = 2, text = "Proportion of Under-Five All-Cause Deaths", outer = TRUE)
    text(x = as.vector(plot)-.3, y = -(.04*ymax), medianGraphs2$variable, adj = c(0,0), cex = .95)
    axis(2, tck=1, ,col.ticks="light gray", labels = FALSE, lty = 2, lwd = 0.75)
    barplot(add = TRUE, medianGraphs2$deathProp,xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = palette[7], las = 1, main = "Median Death Estimates as a Proportion\nof Total Under-Five All-Cause Deaths\n(error bars as Interquartile Range)")
    suppressWarnings(arrows(plot, medianGraphs2$q1/medianGraphs2$whoDeaths, plot, medianGraphs2$q3/medianGraphs2$whoDeaths, angle = 90, code = 3, length = 0.05))

  }
  )

  plotDeathProp2 = function(){
    numReps <- input$N
    results <- generateLHS()
    results <- results[1]; results <- data.frame(results)

    totalDeaths <- results[,40]
    results[, 41] = 1:numReps
    colnames(results) = c(paste(privateSales$Country), "All Countries", "SimulationNumber")

    # reformat data from wide to long format
    meltedOutput <-reshape2::melt(results, id.vars = c("SimulationNumber"))

    # use summary function to generate summary
    meltedSummary <- summarySE(meltedOutput, measurevar = "value", groupvars = c("variable"))

    medianGraphs <- meltedSummary
    medianGraphs2 <- medianGraphs[1:39,]

    # turn names into text
    medianGraphs2$variable <- paste(medianGraphs2$variable)

    # order meltedSummary 2 by country name and load in <5 Deaths
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$variable),]
    medianGraphs2["whoDeaths"] <- WHOdeaths[,3]
    medianGraphs2["whoDeathRate"] <- WHOdeaths[,4]

    # load in <5 malaria deaths
    WHOmalaria <- WHOmalaria[WHOmalaria$Country %in% medianGraphs2$variable,]
    WHOmalaria <- WHOmalaria[order(WHOmalaria$Country),]
    medianGraphs2["malariaDeaths"] <- WHOmalaria$Numeric

    # calculation prop of deaths
    medianGraphs2$whoDeaths <- medianGraphs2$whoDeaths * 1000
    medianGraphs2["malariaProp"] <- medianGraphs2$median / medianGraphs2$malariaDeaths
    medianGraphs2["deathProp"] <- medianGraphs2$median / medianGraphs2$whoDeaths

    #print(medianGraphs2$variable)

    # shorten long country names so that they appear on the graph
    for (i in 1:length(colnames(medianGraphs2))){
      if (medianGraphs2$variable[i] == "Central African Republic") {
        medianGraphs2$variable[i] <- "CAR"
      }
      if (medianGraphs2$variable[i] == "Democratic Republic of the Congo") {
        medianGraphs2$variable[i] <- "DRC"
      }
      if (medianGraphs2$variable[i] == "Cote d'Ivoire"){
        medianGraphs2$variable[i] <- "Coté d'Ivoire"
      }
    }

    #print(medianGraphs2$variable)

    # sort by death prop
    medianGraphs2 <- medianGraphs2[order(medianGraphs2$deathProp, decreasing = TRUE),]

    ymax = max(medianGraphs2$q3/medianGraphs2$whoDeaths)

    #print(ymax)

    par(oma = c(3,1,0,0), srt = -40, xpd = NA)
    plot <-barplot(medianGraphs2$deathProp,xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012), xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = NA, border = NA, las = 1, axes = FALSE)
    mtext(side = 2, text = "Proportion of Under-Five All-Cause Deaths", outer = TRUE)
    text(x = as.vector(plot)-.3, y = -(.02*ymax), medianGraphs2$variable, adj = c(0,0), cex = .85)
    axis(2, tck=1, ,col.ticks="light gray", labels = FALSE, lty = 2, lwd = 0.75)
    barplot(add = TRUE, medianGraphs2$deathProp, xaxs = "i", xlim = xLimsBarPlot(39,padding = 0.012),xlab = NULL, ylab = NULL, ylim = c(0,ymax*1.05), col = palette[7], las = 1, main = "Median Death Estimates as a Proportion\nof Total Under-Five All-Cause Deaths\n(error bars as Interquartile Range)")
    suppressWarnings(arrows(plot, medianGraphs2$q1/medianGraphs2$whoDeaths, plot, medianGraphs2$q3/medianGraphs2$whoDeaths, angle = 90, code = 3, length = 0.05))

  }

  output$deathProp <- renderPlot({
    plotDeathProp()
  }
  )

  # handler for downloading allCauseDeathProp as pdf
  output$downloadDeathProp <- downloadHandler(

    filename = function() {
      results <- generateLHS()
      inputID <- results[5]
      paste(inputID,'_ProportionOfAllCauseDeaths.pdf', sep='')
    }

    ,content = function(file)
    {
      pdf(file, width = 11, height = 8.5)
      plotDeathProp2()
      #plot(1,1)
      #hist(runif(1000, 1, 10))
      dev.off()
    }

    ,contentType = 'application/pdf'
  )

  plotRequestedHist <-reactive({
    country <- input$histogram1
    country2 <- country
    #print(country)

    if (country2 == "Equatorial Guinea") country2 <- "Equatorial.Guinea"
    if (country2 == "Guinea-Bissau") country2 <- "Guinea.Bissau"
    if (country2 == "Central African Republic") country2 <- "Central.African.Republic"
    if (country2 == "Sierra Leone") country2 <- "Sierra.Leone"
    if (country2 == "Burkina Faso") country2 <- "Burkina.Faso"
    if (country2 == "Coté d'Ivoire") country2 <- "Cote.d.Ivoire"
    if (country2 == "Democratic Republic of the Congo") country2 <- "Democratic.Republic.of.the.Congo"
    if (country2 == "All Countries") country2 <- "All.Countries"

    #test <- rnorm(1000, 50, 10)
    #hist(test)

    boxplotOutput <- generateLHS()
    boxplotOutput <- boxplotOutput[1]
    boxplotOutput <- data.frame(boxplotOutput)

    #print(colnames(boxplotOutput))

    LHSsize <- input$N
    #print("LHSsize = ")
    #print(LHSsize)

    for(i in 1:length(colnames(boxplotOutput))){
      boxplotOutput[(LHSsize + 1), i] <- median(boxplotOutput[1:LHSsize, i])
    }

    # order by mean
    boxplotOutput <- boxplotOutput[,order(boxplotOutput[(LHSsize + 1),])]

    # remove mean
    boxplotOutput <- boxplotOutput[1:LHSsize,]

    # print(colnames(boxplotOutput))

    #print(colnames(boxplotOutput))

    firstBox <- boxplotOutput[,country2]
    #firstBox <- as.matrix(firstBox)
    #print(str(firstBox))
    #print("Median:")
    #print(median(boxplotOutput[,1]))
    #print(firstBox)
    maxy <- max(firstBox)

    plot <- hist(firstBox, main = country, xlab = "Estimated Deaths Per Year")
    maxCounts <- max(plot$counts)


    #plot
    if (median(boxplotOutput[,country2]) != 0 | mean(boxplotOutput[,country2]) != 0){
      plot2 <- {
        hist(firstBox, main = country, xlab = "Estimated Deaths Per Year")

        # add red line at the median
        abline(v = median(boxplotOutput[,country2]), lwd = 4, col = "red")

        # add a blue line at the mean
        abline(v = mean(boxplotOutput[,country2]), lwd = 4, col = "blue", lty = 2)

        legend(.85*maxy, maxCounts,c("median","mean"), lty = c(1,2), lwd = c(3,3), seg.len = 3, col = c("red","blue"))
      }
      return(plot2)
    } else{
      plot2 <- {
        hist(firstBox, main = country, xlab = "Estimated Deaths Per Year",xlim = c(-1,10), breaks = c(-.5,.5))

        # add red line at the median
        abline(v = median(boxplotOutput[,country2]), lwd = 4, col = "red")

        # add a blue line at the mean
        abline(v = mean(boxplotOutput[,country2]), lwd = 4, col = "blue", lty = 2)

        legend(6, LHSsize,c("median","mean"), lty = c(1,2), lwd = c(3,3), seg.len = 3, col = c("red","blue"))
      }
      return(plot2)
    }
  }
  )

  plotRequestedHist2 = function()	{
    country <- input$histogram1
    country2 <- country
    #print(country)

    if (country2 == "Equatorial Guinea") country2 <- "Equatorial.Guinea"
    if (country2 == "Guinea-Bissau") country2 <- "Guinea.Bissau"
    if (country2 == "Central African Republic") country2 <- "Central.African.Republic"
    if (country2 == "Sierra Leone") country2 <- "Sierra.Leone"
    if (country2 == "Burkina Faso") country2 <- "Burkina.Faso"
    if (country2 == "Coté d'Ivoire") country2 <- "Cote.d.Ivoire"
    if (country2 == "Democratic Republic of the Congo") country2 <- "Democratic.Republic.of.the.Congo"
    if (country2 == "All Countries") country2 <- "All.Countries"

    #test <- rnorm(1000, 50, 10)
    #hist(test)

    boxplotOutput <- generateLHS()
    boxplotOutput <- boxplotOutput[1]
    boxplotOutput <- data.frame(boxplotOutput)

    #print(colnames(boxplotOutput))

    LHSsize <- input$N
    #print("LHSsize = ")
    #print(LHSsize)

    for(i in 1:length(colnames(boxplotOutput))){
      boxplotOutput[(LHSsize + 1), i] <- median(boxplotOutput[1:LHSsize, i])
    }

    # order by mean
    boxplotOutput <- boxplotOutput[,order(boxplotOutput[(LHSsize + 1),])]

    # remove mean
    boxplotOutput <- boxplotOutput[1:LHSsize,]

    # print(colnames(boxplotOutput))

    #print(colnames(boxplotOutput))

    firstBox <- boxplotOutput[,country2]
    #firstBox <- as.matrix(firstBox)
    #print(str(firstBox))
    #print("Median:")
    #print(median(boxplotOutput[,1]))
    #print(firstBox)
    maxy <- max(firstBox)


    if (median(boxplotOutput[,country2]) != 0 | mean(boxplotOutput[,country2]) != 0){
      plot <- hist(firstBox, main = country, xlab = "Estimated Deaths Per Year")
      maxCounts <- max(plot$counts)

      # add red line at the median
      abline(v = median(boxplotOutput[,country2]), lwd = 4, col = "red")

      # add a blue line at the mean
      abline(v = mean(boxplotOutput[,country2]), lwd = 4, col = "blue", lty = 2)

      legend(.85*maxy, maxCounts,c("median","mean"), lty = c(1,2), lwd = c(3,3), seg.len = 3, col = c("red","blue"))

    } else{
      plot <- hist(firstBox, main = country, xlab = "Estimated Deaths Per Year",xlim = c(-1,10), breaks = c(-.5,.5))

      # add red line at the median
      abline(v = median(boxplotOutput[,country2]), lwd = 4, col = "red")

      # add a blue line at the mean
      abline(v = mean(boxplotOutput[,country2]), lwd = 4, col = "blue", lty = 2)

      legend(6, LHSsize,c("median","mean"), lty = c(1,2), lwd = c(3,3), seg.len = 3, col = c("red","blue"))
    }
  }

  # handler for downloaded requested Histogram as pdf
  output$downloadHist <- downloadHandler(

    filename = function() {
      results <- generateLHS()
      inputID <- results[5]
      #print(inputID)
      paste(inputID,'_Histogram_',input$histogram1,'.pdf', sep='')
      #name <- "HistogramTest.png"
      #name
    }

    ,content = function(file)
    {
      pdf(file, width = 11, height = 8.5)
      plotRequestedHist2()
      #plot(1,1)
      #hist(runif(1000, 1, 10))
      dev.off()
    }

    ,contentType = 'application/pdf'
  )

  # test = function(){
  # nameTEST <- c("HistogramTest1.png")
  # #print(nameTEST)

  # png(nameTEST)

  # plotRequestedHist()
  # dev.off()
  # #print("finished")
  # }

  # return the requested histogram
  output$requestedHist <- renderPlot({
    plotRequestedHist()
  }
  )

  generatePRCC <- reactive({
    out <- generateLHS()
    sameFlag <- out[6]
    #print(sameFlag)
    if(sameFlag == TRUE){
      stop("Partial Rank Correlation Coefficients cannot be calculated when at least one of the input parameters does not have a range (i.e. the min and max slider for that input parameter are set to the same value).")
    }
    results <- out[1]
    results <- data.frame(results)

    deathRate <- out[2]
    deathRate <- data.frame(deathRate)

    sales <- out[3]
    sales <- data.frame(sales)

    fakePercent <-out[4]
    fakePercent <- data.frame(fakePercent)

    totalDeaths = results[,40]

    names2 <- names[-40]
    names2 <- paste(names2)

    #print(class(names2))

    # shorten long country names so that they appear on the graph
    for (i in 1:length(names2)){
      if (names2[i] == "Central African Republic") {
        names2[i] <- "CAR"
      }
      if (names2[i] == "Democratic Republic of the Congo") {
        names2[i] <- "DRC"
      }
      if (names2[i] == "Cote d'Ivoire"){
        names2[i] <- "Coté d'Ivoire"
      }
    }

    salesLabel <- replicate(39, ": Antimalarial Sales")
    salesLabel <- paste(paste(names2), salesLabel, sep = "")

    fakeLabel<- replicate(39, ": Prevalence of PQ Antimalarials")
    fakeLabel <- paste(paste(names2), fakeLabel, sep = "")

    x.base <- data.frame(cbind(sales, fakePercent, deathRate))
    x = data.frame(x.base, totalDeaths)

    colnames(x) <- c(salesLabel, fakeLabel, "Case Fatality Rate", "TotalDeathsOutput")
    totDeathSens <- counterfeitPRCC(x, sort.results = TRUE, sort.abs = TRUE)
    totDeathSens <- totDeathSens[-2]

    totDeathSens <- data.frame(rownames(totDeathSens), totDeathSens[,1], totDeathSens[,2], stringsAsFactors = FALSE)
    colnames(totDeathSens) <- c("Input Parameter","PRCC", "P-value")
    totDeathSens <- totDeathSens[order(totDeathSens$PRCC, decreasing = TRUE),]
    rownames(totDeathSens) <- NULL

    #print(class(totDeathSens$PRCC))

    #print(totDeathSens)
    # return PRCC table
    totDeathSens
  }
  )

  # generate an HTML table view of the PRCCs
  output$PRCC <- renderTable({
    data.frame(generatePRCC())
  }
  )

  # function to allow for .csv download of PRCCs
  output$downloadPRCC <- downloadHandler(
    filename = function() {
      results <- generateLHS()
      inputID <- results[5]
      #print(inputID)
      paste(inputID,'_SensitivityAnalysis.csv', sep='')
    }
    ,content = function(file){
      out <- generatePRCC()
      #print(class(out[,1]))

      for(i in 1:length(out[,1])){
        if(out[i,1] == "Coté d'Ivoire: Antimalarial Sales"){
          out[i,1] <- "Cote d'Ivoire: Antimalarial Sales"
        }
        if(out[i,1] == "Coté d'Ivoire: Prevalence of PQ Antimalarials"){
          out[i,1] <- "Cote d'Ivoire: Prevalence of PQ Antimalarials"
        }
      }

      write.csv(out, file)
    }
  )

  generateInputs <-reactive({
    countrySpec <- input$countrySpecific
    manuscript <- input$manuscript
    names2 <- names[-40]
    names2 <- paste(names2)

    #print(class(names2))

    # shorten long country names so that they appear on the graph
    for (i in 1:length(names2)){
      if (names2[i] == "Central African Republic") {
        names2[i] <- "CAR"
      }
      if (names2[i] == "Democratic Republic of the Congo") {
        names2[i] <- "DRC"
      }
      #if (names2[i] == "Cote d'Ivoire"){
      #	names2[i] <- "Coté d'Ivoire"
      #}
    }

    displaySales <- data.frame(displaySales)
    rownames(displaySales) <- NULL

    #
    salesLabel1 <- replicate(39, "Antimalarial Sales")
    salesLabel2 <- paste(names2)
    salesMin <- replicate(39, "NA")
    salesMax <- replicate(39, "NA")
    salesDist <- replicate(39, "Normal")

    #
    fakeLabel1<- replicate(39, "Prevalence of PQ Antimalarials")
    fakeLabel2 <- paste(names2)
    fakeDist <- replicate(39, "Uniform")

    salesMean <- replicate(39, 0)
    salesStd <- replicate(39,0)

    fakeMean <- replicate(39,0)
    fakeStd <- replicate(39,"NA")

    if(countrySpec) {
      # make dependent on Update button
      input$goButton

      fakeMin <- isolate(c(input$Prev1[1],input$Prev2[1],input$Prev3[1],input$Prev4[1],input$Prev5[1],input$Prev6[1],input$Prev7[1],input$Prev8[1],input$Prev9[1],input$Prev10[1],input$Prev11[1],input$Prev12[1],input$Prev13[1],input$Prev14[1],input$Prev15[1],input$Prev16[1],input$Prev17[1],input$Prev18[1],input$Prev19[1],input$Prev20[1],input$Prev21[1],input$Prev22[1],input$Prev23[1],input$Prev24[1],input$Prev25[1],input$Prev26[1],input$Prev27[1],input$Prev28[1],input$Prev29[1],input$Prev30[1],input$Prev31[1],input$Prev32[1],input$Prev33[1],input$Prev34[1],input$Prev35[1],input$Prev36[1],input$Prev37[1],input$Prev38[1],input$Prev39[1]))

      fakeMax <- isolate(c(input$Prev1[2],input$Prev2[2],input$Prev3[2],input$Prev4[2],input$Prev5[2],input$Prev6[2],input$Prev7[2],input$Prev8[2],input$Prev9[2],input$Prev10[2],input$Prev11[2],input$Prev12[2],input$Prev13[2],input$Prev14[2],input$Prev15[2],input$Prev16[2],input$Prev17[2],input$Prev18[2],input$Prev19[2],input$Prev20[2],input$Prev21[2],input$Prev22[2],input$Prev23[2],input$Prev24[2],input$Prev25[2],input$Prev26[2],input$Prev27[2],input$Prev28[2],input$Prev29[2],input$Prev30[2],input$Prev31[2],input$Prev32[2],input$Prev33[2],input$Prev34[2],input$Prev35[2],input$Prev36[2],input$Prev37[2],input$Prev38[2],input$Prev39[2]))

    } else{
      fakeMin <- replicate(39,input$Prev[1])
      fakeMax <- replicate(39,input$Prev[2])
    }

    j = 1
    surveyCountries = c(29,33,32,35,25,39,15,34,37)
    for(i in 1:39){
      # Sales Inputs
      salesMean[i] <- displaySales$Mean[i]
      salesStd[i] <- displaySales$Standard_Deviation[i]

      # Prevalence Inputs
      if(manuscript){
        if(i %in% surveyCountries){
          fakeDist[i] <- "Normal"
          fakeMean[i] <- antimalarialQualitySurvey$Mean[j]
          fakeStd[i] <- antimalarialQualitySurvey$SE[j]
          fakeMin[i] <- "NA"
          fakeMax[i] <- "NA"
          j = j + 1
        } else {
          fakeMean[i] <- (0.4+0)/2
          fakeStd[i] <- "NA"
        }
      } else {
        fakeMean[i] <- (fakeMax[i]+fakeMin[i])/2
        fakeStd[i] <- "NA"
      }
    }

    # create and sort sales inputs by name
    cols <- c("Country","InputParameter","Distribution Shape", "Min", "Mean","Max", "StdDev")
    salesINPUTS <- data.frame(salesLabel2,salesLabel1,salesDist, salesMin, salesMean, salesMax, salesStd, stringsAsFactors = FALSE)
    colnames(salesINPUTS) <- cols
    #print(salesINPUTS)
    salesINPUTS <- salesINPUTS[order(salesINPUTS$Country),]
    cols <- c("Country","Input Parameter","Distribution Shape", "Min", "Mean","Max", "StdDev")
    colnames(salesINPUTS) <- cols

    # create and sort PQ prevalence inputs by name
    cols <- c("Country","InputParameter","Distribution Shape", "Min", "Mean","Max", "StdDev")
    fakeINPUTS <- data.frame(fakeLabel2,fakeLabel1, fakeDist, fakeMin, fakeMean, fakeMax, fakeStd, stringsAsFactors = FALSE)
    colnames(fakeINPUTS) <- cols
    #print(fakeINPUTS)
    fakeINPUTS <- fakeINPUTS[order(fakeINPUTS$Country),]
    cols <- c("Country","Input Parameter","Distribution Shape", "Min", "Mean","Max", "StdDev")
    colnames(fakeINPUTS) <- cols

    #case fatality rate inputs
    cfrINPUTS <- list("All Countries","Case Fatality Rate","Uniform", input$CFR[1], (input$CFR[2]+input$CFR[1])/2,
                      input$CFR[2], "NA")

    masterINPUTS <- rbind(cfrINPUTS, fakeINPUTS, salesINPUTS)

    masterINPUTS$Min <- as.character(masterINPUTS$Min)
    masterINPUTS$Mean <- as.character(masterINPUTS$Mean)
    masterINPUTS$Max <- as.character(masterINPUTS$Max)
    masterINPUTS$StdDev <- as.character(masterINPUTS$StdDev)

    masterINPUTS
  }
  )

  output$InputParameters <- renderTable({
    inputs <- data.frame(generateInputs())
    rownames(inputs) <- NULL
    inputs
  })

  output$downloadInputs <- downloadHandler(
    filename = function() {
      results <- generateLHS()
      inputID <- results[5]
      #print(inputID)
      paste(inputID,'_InputParamters.csv', sep='')
    }
    ,content = function(file){
      out <- generateInputs()
      rownames(out) <- NULL
      #print(out)
      for(i in 1:length(out[,1])){
        if(out[i,1] == "Coté d'Ivoire: Antimalarial Sales") out[i,1] <- "Cote d'Ivoire: Antimalarial Sales"
        if(out[i,1] == "Coté d'Ivoire: Prevalence of PQ Antimalarials") out[i,1] <- "Cote d'Ivoire: Prevalence of PQ Antimalarials"
      }
      write.csv(out, file)
    }
  )


  # generate an HTML Table view of the Cohen Antimalarial Sales
  # input data
  output$sales <- renderTable({
    #print("I AM HERE Third")
    data.frame(displaySales)
    displaySales <- displaySales[order(displaySales$Mean, decreasing = TRUE),]
    rownames(displaySales) <- NULL
    displaySales
  }
  )

}
)
