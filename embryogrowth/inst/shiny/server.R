library(shiny)
# runApp(".", launch.browser = TRUE)

library(mapdata)
library(embryogrowth)

# Define server logic required to plot various variables against mpg

shinyServer(function(input, output, session) {
  
  observe({
    # file1_NestingArea <- STSRE_NestingArea
    
    output$countrySelect=renderUI({
      countryChoices=subset(STSRE_NestingArea, Sp==input$species)$Country
      selectInput("country", "Country:", choices=countryChoices)
    })
    
    output$MAP=renderPlot({
      NestingArea=c(x=as.numeric(gsub(",", ".", input$longitude)),
                    y=as.numeric(gsub(",", ".", input$latitude)))
      
      zoom=as.numeric(input$zoom)*10
      
      a=map('worldHires', plot=FALSE)
      
      distance=(a$x-NestingArea["x"])^2+(a$y-NestingArea["y"])^2
      
      d=which.min(distance)
      
      NestingArea.coord=c(x=a$x[d], y=a$y[d])
      
      Map=map('worldHires', names=TRUE, xlim=c(NestingArea.coord["x"], NestingArea.coord["x"]), 
              ylim=c(NestingArea.coord["y"], NestingArea.coord["y"]), plot=FALSE)
      
      map('worldHires', xlim=c(NestingArea["x"]-zoom, NestingArea["x"]+zoom), 
          lwd=1, fill=TRUE, col="green", asp=1)
      
      points(NestingArea["x"], NestingArea["y"], pch="o", col="red")
      
      output$NestingAreaMap=renderText(Map)
    })
    
    
    output$rmuSelect=renderUI({
       rmuChoices=subset(STSRE_NestingArea, Sp==input$species & Country==input$country)$RMU
      if (length(rmuChoices)==0) rmuChoices <- "Not possible"
      selectInput("rmu", "RMU:", choices=rmuChoices)
    })
    
    output$resultPlot=renderPlot({
      
      # input$file1 will be NULL initially. After the user selects and uploads a 
      # file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
      # columns. The 'datapath' column will contain the local filenames where the 
      # data can be found.
      
      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      
      if (length(input$rmu)==0) {
        output$resultsInfo <- renderText("Choose first the RMU of your population")
        return(NULL)
      }
      
      modelgrowth <- NULL
      modelgrowth.se <- NULL
      newp <- NULL
      
      if (input$species=="Cc") {
        modelgrowth <- resultNest_4p
        newp <- GenerateAnchor(nests=resultNest_4p, number.anchors=7)
        modelgrowth.se <- result_mcmc_newp$TimeSeriesSE
      }
      
      if (is.null(modelgrowth)) {
        output$resultsInfo <- renderText("Embryonic growth rate is not available for this species")
        return(NULL)
      }
      
      # dest <- "/Users/marc/Documents/Espace_de_travail_R/shiny/STSRE_MSM/Essai.csv"
      # table <- read.csv(dest, header=TRUE, sep=";", stringsAsFactors=FALSE)
      table <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
      if (dim(table)[2]!=2) {
        output$resultsInfo <- renderText("Check the data; 2 columns are required")
        return(NULL)
      }
      
      time <- strptime(table[,1], format=input$dateformat)
      # time <- strptime(table[,1], format='%d/%m/%Y %H:%M:%S')
      time.min <- floor(as.numeric(time-time[1])/60)
      df.nf <- data.frame(Time=time.min, Temperature=table[,2])
      formated <- FormatNests(df.nf)
      out <- info.nests(modelgrowth, parameters=newp, SE=modelgrowth.se, replicate.CI = 10,
                        temperatures=formated, stopattest=TRUE, progress=FALSE)
      output$resultsInfo <- renderText(paste0("species: ", switch(input$species, 
                                                                 Cc="Caretta caretta", 
                                                                 Dc="Dermochelys coriacea",
                                                                 Cm="Chelonia mydas",
                                                                 Ei="Eretmochelys imbricata",
                                                                 Lo="Lepidochelys olivacea",
                                                                 Lk="Lepidochelys kempii",
                                                                 Nd="Natator depressus"), 
                                              "\nRMU: ", input$rmu,
                                              "\nTemperature average: ", mean(table[,2]),
                                              "\nMaximal time in minutes: ", tail(time.min, n=1),
                                              "\nCTE using weighted growth: ", as.numeric(out$summary[[1]]["weighted.temperature.mean"]),
                                                 "SE ", as.numeric(out$summary[[1]]["weighted.temperature.SE"]),
                                              "\nIncubation duration in days: ", as.numeric(out$summary[[1]]["incubation.length.mean"]),
                                                 "SE ", as.numeric(out$summary[[1]]["incubation.length.SE"])))
        plotR(modelgrowth, SE=modelgrowth.se, parameters=newp, replicate.CI=10)
    })
    
    output$TSDPlot=renderPlot({
      csttemp <- data.frame(Temp=numeric(), Total=integer(), Males=integer(), Females=integer())
      for (ref in levels(as.factor(STSRE_TSD$Ref)))      
        if (input[[paste0("c", ref)]]) csttemp <- rbind(csttemp, subset(STSRE_TSD, Ref==ref))
      
      if (dim(csttemp)[1]!=0) {
        rtsd <- tsd(males=csttemp$Males, females=csttemp$Females, temperatures=csttemp$Temp, print=FALSE)
        output$TSDInfo <- renderText(paste("The Pivotal temperature is", sprintf("%.3f",as.numeric(rtsd$par["P"])),
"SE", sprintf("%.3f",as.numeric(rtsd$SE["P"])), "\nThe Transitional range of temperatures l=5% is", sprintf("%.3f",as.numeric(rtsd$TRT)),
"SE", sprintf("%.3f",as.numeric(rtsd$SE_TRT))
))
        } else {
        output$TSDInfo <- renderText('No data have been selected')
      }
    })
        
    if (!is.null(input$click)) {
      updateNumericInput(session, "latitude", value = gsub(",", ".", as.character(input$click$y)))
      updateNumericInput(session, "longitude", value = gsub(",", ".", as.character(input$click$x)))
    }
  })
  
})
