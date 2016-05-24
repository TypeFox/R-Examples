AdaptGauss = function(Data,Means=NaN,SDs=NaN,Weights=NaN,ParetoRadius=NaN,LB=NaN,HB=NaN,ListOfAdaptGauss,fast=TRUE){
#    Out=AdaptGauss2(Data,Means,SDs,Weights,ParetoRadius,LB,HB);
#    adapt interactively a Gaussians Mixture Model GMM to the empirical PDF of the data such that N(M,S)*W is a  model for Data
#   
#    INPUT
#    Data(1:n)                      Vector of Data, may contain NaN
#   
#    OPTIONAL
#    Means(1:c)                     The means of the distribution;             default: nanmean(Data)
#    SDs(1:c)                       The StandardDeviatons   of the distributions            default: nanstd(Data)
#    Weights(1:c)                   The weights of the distributions           default: 1
#    ParetoRadius                   the ParetoRadius for PDE on Data,
#                                   It is calculated, if not given
#    LB,HB                          limits where the adaptation is done, default: [min(Data) max(Data)]
#    ListOfAdaptGauss               if given in return Format of AdaptGauss, the modell can be edited
#
#    OUTPUT
#    Out$Means(1:L)                 The adapted means of the distributions
#    Out$SDs(1:L)                   The adapted sdev  of the distributions              
#    Out$Weights(1:L)               The adapted weights of the distributions
#    Out$ParetoRadius               Pareto Radius used for empirical PDE 
#    Out$RMS                        Root Mean Square Distance between empirical PDE and pdf(GMM)
#                                   RMS == sqrt(sum(empirical PDE - pdf(GMM))^2)
#   Out$BayesBoundaries             Bayes Boundaries between gaussians
# author: Onno-Hansen Goos
# 1.Editor: MT 08/2015
 Data; #Bricht bei nicht existierendem Bezeichner ab
  ## Starte Shiny
  #library(shiny)

if(!missing(ListOfAdaptGauss)){
  if(is.list(ListOfAdaptGauss)){
    Means<-ListOfAdaptGauss$Means
    if(is.null(Means)) Means=NaN
    SDs<-ListOfAdaptGauss$SDs
    if(is.null(SDs)) SDs=NaN
    Weights<-ListOfAdaptGauss$Weights
    if(is.null(Weights)) Weights=NaN
  }
}

  if (!(hasArg(Data))){
    stop("No Data Input.")
   }  

  if(!missing(Data)){
      if(!is.vector(Data)){
        stop("Data has to be a vector, maybe use as.vector(Data)")
      }
     if(!is.numeric(Data)){
         stop("Data has to be numeric, maybe use as.numeric(Data)")     
     }
  }  
  
  # rsampleAdaptGauss
  ############################################################
  `rsampleAdaptGauss` <- function(k,n,uniq=TRUE,exact=TRUE){
    
    index = ceiling(runif(k)*n) # Calculate k of n values.
    if(uniq==TRUE){# if unique, avoide duplicates.
      if(k>n&exact){ # Don't run in infinite loop.
        print("The first parameter has to be <= the second.")
        print("You choose uniq = TRUE. This would cause an infinite loop.")
        print("Abort function.")
      }
      else{
        index = unique(index)
        while(exact&(length(index)<k)){
          index = c(index,ceiling(runif(k-length(index))*n))
          index = unique(index)
        }
        return(index)
      }
    }
    else{
      return(index)
    }
  }
  ############################################################
  
  # getOptGauss
  #########################################################################
  `getOptGauss` <- function(Data, Kernels, ParetoDensity,fast){ 
    # Teste RMS fuer einen Gauss
    Mean1 <- mean(Data)
    Deviation1 <- sd(Data)
    Weight1 <- 1
    Var=EMGauss(Data,fast=fast)
    Mean1 <- Var$Means
    Deviation1 <- Var$SDs
    Weight1 <- Var$Weights
    Fi <- dnorm(Kernels,Mean1,Deviation1)*Weight1
    RMS1 <- sqrt(sum((Fi-ParetoDensity)^2))
    
    # Teste RMS fuer 3 Gauss
    Means2 <- c(0,0,0)
    Deviations2 <- c(0,0,0)
    Weights2 <- c(0,0,0)
    Valskmeans <- kmeans(Data,3,iter.max=100)
    KValues <- Valskmeans$cluster
    #print(KValues2)
    for (i in 1:3){
      Means2[i] <- mean(Data[KValues==i])
      Deviations2[i] <- sd(Data[KValues==i])
      Weights2[i] <- sum(KValues==i)
      if (is.na(Deviations2[i])) {Deviations2[i] <- 0}
    }
    Weights2 <- Weights2/length(KValues)
    Var=EMGauss(Data,Means2,Deviations2,Weights2,10,fast=fast)
    Means2 <- Var$Means
    Deviations2 <- Var$SDs
    Weights2 <- Var$Weights
    Fi <- 0
    for (i in 1:3){
      Fi <- Fi+dnorm(Kernels,Means2[i],Deviations2[i])*Weights2[i]
    }
    RMS2 <- sqrt(sum((Fi-ParetoDensity)^2))
    
    # ueberpruefe ob RMS1( 1 Gauss) oder RMS2 (3 Gauss ) kleiner ist. Speichere zugehoerige means, deviations und weights
    SSE1=RMS1^2*log(3)
    SSE2=RMS2^2*log(3*3)
    if (SSE1<SSE2){
      means <- Mean1
      deviations <- Deviation1
      weights <- Weight1
    } else {
      means <- Means2
      deviations <- Deviations2
      weights <- Weights2
    }
    # Ordne gaussians nach mean
    order <- order(means)
    means <- means[order]
    deviations <- deviations[order]
    weights <- weights[order]
    out=list(means=means,deviations=deviations,weights=weights)
    
    
    return(out)
  }
  #########################################################################
  


  outputApp=runApp(list(
    ## ui.R --> stellt oberfl?che her
    ui = fluidPage(
      #Ueberschrift-Panel
      #headerPanel("Adapt Gauss"), # oh! die Kommas sind wichtig! ^^
    fluidRow(
      column(3,
               fluidRow(
                 wellPanel(# Legt aktuellen Gauss fest (welcher gerade bearbeitet werden kann)
                   tags$div(class = "row",
                   uiOutput("sliderCurrGauss"),
                   actionButton("AddGaussButton", "Add", icon = icon("plus-square")),
                   actionButton("DeleteGaussButton", "Delete Current",icon = icon("bitbucket"))#,
                   )
                 )
              )        
      ),
      column(2,offset = 1,
                            fluidRow(
                            h6("Expectation Maximation Algorithm with Iterations:"),
                            numericInput("numIterations", 
                                          label ="", value = 10),
                             actionButton("expMax", "", icon = icon("calculator")),
                            actionButton("restore", " ", icon = icon("undo"))
                          )       
      ), 
      column(3,
                          #wellPanel(# weights
                            h6("Weights"),
                            actionButton("normAll", "Normalize All", icon = icon("balance-scale")),
                            actionButton("normOth", "Normalize Other", icon = icon("angle-double-down")),
                            h6("Options"),
                            checkboxInput("showComponents", "Show Components", TRUE),
                            checkboxInput("showBayes", "Show Bayes Boundaries", FALSE)
                          #)
      ),
      column(3, 
                   wellPanel(
                    h6("General"),
                    actionButton("PlotFig", "Plot Figure", icon = icon("area-chart")),
                    actionButton("RestoreBestRMS", "Restore Best RMS", icon = icon("backward")),           
                    h6("AdaptGauss:"),
                    actionButton("CloseButton", "Close", icon = icon("close"))
                   )
      ),
      fluidRow(
                   column(12,offset=1,
                      plotOutput('PDE',width = 700, height = 450)
                   )
      ),
      fluidRow(
                   column(2,
                    wellPanel(h5(textOutput('gaussianNumber'))),
                    uiOutput("minMean"),
                    uiOutput("minSdev"),
                    uiOutput("minWeight")
      ),
                   column(8,
                     uiOutput("sliderM"),
                     uiOutput("sliderS"),
                     uiOutput("sliderW")
                   ),
                   column(2,
                      wellPanel(h5(textOutput('ScreenMessage'))),
                      uiOutput("maxMean"),
                      uiOutput("maxSdev"),
                      uiOutput("maxWeight")
                   )
      )
    )),
    server = function(input, output, session){
      ## server.R --> Ab hier keine Oberfl?che mehr, ist aber mit ui.R (Oberfl?che) verkn?pft
      
      
      # Default values for Data if no input

      data <- Data
      
      #?bertrage Input Werte (k?nnen sonst von Shiny nicht korrekt verwendet werden)
      
      GM <- Means
      GS <- SDs
      GW <- Weights
      ParetoRadius <- ParetoRadius
      LB <- LB
      HB <- HB
      nSignif <- 4 # Auf wie viele "echte" Stellen soll nach dem Komma gerundet werden (zB. RMS)

      numGaussSave <- NULL
      MLimit <- NULL
      GSSave <- NULL
      GSBestRMS <- NULL
      GMSave <- NULL
      GMBestRMS <- NULL
      GWSave <- NULL
      GWBestRMS <- NULL

      
      # Index f?r Befehle (siehe interactive value befehl)
      iBefehl <- 0

      # Lada Daten, definiere Variablen
      #observe({ # Wird ein mal am Anfang ausgef?hrt. Erzeugt means, stds usw.
      print("Start Session")
      #print("Load Data")        
      if(length(data)==0)
        return("Error: Data could not be loaded");
      
      # Wie viele Gauss sollen getestet werden? noch nicht implentiert!
      suggestedMaxGauss <- 5
      
      # Eliminate NaNs in Data 
      dataNew <- 0
      j <- 1
      for (i in 1:length(data)){
        if (!is.nan(data[i]) && !is.na(data[i])){
          dataNew[j] <- data[i]
          j <- j+1
        }
      } 
      data <- dataNew
      # Setzte LB (low Boundary) und HB (high boundary) Falls nicht im Input
      if(is.nan(LB)) LB <- min(data)
      if(is.nan(HB)) HB <- max(data)
      # Eliminate Values below LB and above HB in Data
      dataNew <- 0
      j <- 1
      for (i in 1:length(data)){
        if (data[i]>=LB && data[i]<=HB){
          dataNew[j] <- data[i]
          j <- j+1
        }
      } 
      data <- dataNew
      # Reduce Data to 10000 Elements, if larger than 10000 (?bersch?ssige Datenpunkte werden randomisiert entfernt)
      if (length(data)>10000){
        data <- data[rsampleAdaptGauss(10000,length(data))];
        print("Reducing to 10000 datapoints");
      } 
      # Bestimme Pareto Density inkl. Kernels
      if (is.nan(ParetoRadius)){ 
        ParetoRadius<-ParetoRadius(data) 
        nRow=length(data)
        #MT: Halte ich nicht fuer plausibel
        #if (nRow>1024){
       #   ParetoRadius = ParetoRadius * 4 /(nRow^0.2);
        #}
      }
        ParetoDensityEstimationVar <- ParetoDensityEstimation(data,paretoRadius=ParetoRadius);
      ParetoDensity <- ParetoDensityEstimationVar$paretoDensity;
      Kernels <- ParetoDensityEstimationVar$kernels; 
      # Setze Werte f?r Means, deviations und weights, falls nicht ?bergeben
      if (is.nan(sum(GM)) || is.nan(sum(GS)) || is.nan(sum(GW)) || length(GM)!=length(GS) || length(GM)!=length(GW) ){
        vars=getOptGauss(Data=data, Kernels, ParetoDensity,fast=fast)
        GM <- vars$means
        GS <- vars$deviations
        GW <- vars$weights
      }
      BB <- NaN # Bayes Boundaries (wird berechnet, bevor es das erste mal ausgegeben wird)
      #GM <<- means          #Mean der Gaussians
      #GS <<- deviations          #StdDev der Gaussians
      #GW <<- weights        #Weight der Gaussians 
      numGauss <- length(GM)#Anzahl der Gaussians
      numIterations <- 10  # Default value for number of Iterations (in EMGauss)
      RMS <- 99 # Root Mean Sqare (wird berechnet, bevor es das erste mal ausgegeben wird)
      
      meanRMS0 <- mean(data)
      DeviationRMS0 <- sd(data)
      Fi <- dnorm(Kernels,meanRMS0,DeviationRMS0)
      RMS0 <- sqrt(sum((Fi-ParetoDensity)^2))
      
      currGauss <- 1 # Current Gauss (der Gauss, welcher gerade bearbeitet werden kann)
      # Speicher Werte f?r "restore Best RMS
      GMBestRMS <- GM
      GSBestRMS <- GS
      GWBestRMS <- GW
      BestRMS <- RMS
      numGaussBestRMS <- numGauss
      currGaussBestRMS <- currGauss
      
      xlimit <- c(min(Kernels), max(Kernels));
      ylimit <- c(0,max(ParetoDensity*1.2))
      numKernels <- length(Kernels)
      MDefault <- mean(data) # Default Werte f?r neu erzeugte Gauss
      SDefault <- sd(data)
      WDefault <- 0.5
      MLimit <- xlimit # Legt bereiche fest, in welchen sich die Werte f?r die Gaussians befinden k?nnen
      SLimit <- c(0,max(SDefault*2,max(GS)))
      WLimit <- c(0,1)
        #print("[DONE]");
      #})

      
      ## Output Part: lege slider, Felder und Buttons fest (sind im ui.R part eingebunden)
      # Text Output
      output$Header <- renderText({'PDE Plot of uploaded Data'})
      output$gaussianNumber <- renderText({
        befehl$updateCurrGauss
        paste0("GaussianNo.",currGauss)
      })
      # Slider Output
      output$sliderCurrGauss <- renderUI({
        if (numGauss > 1){
          sliderInput('currGauss', 'Gaussian No.',ticks=FALSE,width='50%', sep="" ,min = 1, max = numGauss, value = currGauss, step= 1)
        }
      })
      output$sliderM <- renderUI({
        befehl$drawSliderMSW
        sliderInput("M", h5("Mean"), width='150%', min=MLimit[1], max=MLimit[2], value = GM[currGauss], step= (MLimit[2]-MLimit[1])*0.001)
      })  
      output$sliderS <- renderUI({
        befehl$drawSliderMSW
        sliderInput("S", h5("SD"), width='150%', min=SLimit[1], max=SLimit[2], value = GS[currGauss], step=(SLimit[2]-SLimit[1])*0.001)
      })  
      output$sliderW <- renderUI({
        befehl$drawSliderMSW
        sliderInput("W", h5("Weight"), width='150%',min=WLimit[1], max=WLimit[2], value = GW[currGauss], step=(WLimit[2]-WLimit[1])*0.001)
      })       
      output$minMean <- renderUI({
        befehl$updateMinMax
        #print("Erstelle Min/Max f?r M/S/W")
        numericInput("minMean", 
                     label = h6("Min Mean"), 
                     value = signif(MLimit[1],digits=nSignif)
        )
      }) 
      output$maxMean <- renderUI({
        befehl$updateMinMax
        numericInput("maxMean", 
                     label = h6("Max Mean"), 
                     value = signif(MLimit[2],digits=nSignif)
        )
      })
      output$minSdev <- renderUI({
        befehl$updateMinMax
        numericInput("minSdev", 
                     label = h6("Min SD"), 
                     value = signif(SLimit[1],digits=nSignif)
        )
      }) 
      output$maxSdev <- renderUI({
        befehl$updateMinMax
        numericInput("maxSdev", 
                     label = h6("Max SD"), 
                     value = signif(SLimit[2],digits=nSignif)
        )
      }) 
      output$minWeight <- renderUI({
        befehl$updateMinMax
        numericInput("minWeight", 
                     label = h6("Min Weight"), 
                     value = signif(WLimit[1],digits=nSignif)
        )
      }) 
      output$maxWeight <- renderUI({
        befehl$updateMinMax
        numericInput("maxWeight", 
                     label = h6("Max Weight"), 
                     value = signif(WLimit[2],digits=nSignif)
        )
      }) 
      
      
      ## Interaktive Variablen: befehl$... fungiert als Funktionsaufruf.
      befehl <- reactiveValues(plot = 0, updateSlider=0, updateSliderCurrGauss=0, drawSliderMSW=0, updateRMS=0, updateMinMax=0) #Signal zum update des Plots / der Slider
      
      
      ## Observe Part: warte auf Input     
      observe({ # Grenzen f?r die Werte von means, sdevs und weights
        if (!is.null(input$minMean)){
          if (  is.numeric(input$minMean)  &&  input$minMean<min(GM)  )  MLimit[1] <<- input$minMean
          if (  is.numeric(input$maxMean)  &&  input$maxMean>max(GM)  )  MLimit[2] <<- input$maxMean
          if (  is.numeric(input$minSdev)  &&  input$minSdev<min(GS)  )  SLimit[1] <<- input$minSdev
          if (  is.numeric(input$maxSdev)  &&  input$maxSdev>max(GS)  )  SLimit[2] <<- input$maxSdev
          if (  is.numeric(input$minWeight)  &&  input$minWeight<min(GW)  )  WLimit[1] <<- input$minWeight
          if (  is.numeric(input$maxWeight)  &&  input$maxWeight>max(GW)  )  WLimit[2] <<- input$maxWeight
          iBefehl <<- iBefehl+1
          befehl$drawSliderMSW <- iBefehl
        }
      })
      
      # Normalisiert alle Gaussians
      observe({
        #print("Normalize All")
        if (input$normAll>0){
          sumWeight <- sum(GW)
          for (i in 1:numGauss){
            GW[i] <<- GW[i]/sumWeight
          }
          iBefehl <<- iBefehl+1
          befehl$plot <- iBefehl
          befehl$updateSlider <- iBefehl
        }
      })
      # Normalisiert andere Gaussians (alle ausser dem aktuellen)
      observe({
        #print("Normalize Other")
        if (input$normOth>0){
          sumWeight <- sum(GW)
          sumWeightOther <- sumWeight-GW[currGauss]
          zielSum <- 1-GW[currGauss]
          for (i in 1:numGauss){
            if (i!=currGauss){
              GW[i] <<- GW[i]/sumWeightOther*zielSum
            }
          }
          iBefehl <<- iBefehl+1
          befehl$plot <- iBefehl
          befehl$updateSlider <- iBefehl
        }
      })
      # wechsel den aktuellen Gauss
      observe({
        #print("Update currGauss")
        if (!is.null(input$currGauss)){
          currGauss <<- input$currGauss
          iBefehl <<- iBefehl+1
          befehl$updateCurrGauss <- iBefehl
        }
      })
      #Erneuert Slider f?r M/S/W, (zB. wenn der Gauss gewechselt wurde)
      observe({
        befehl$updateSlider
        befehl$updateCurrGauss
        #print("Update Slieder M/S/W")
        updateSliderInput(session, "M", value = GM[currGauss])
        updateSliderInput(session, "S", value = GS[currGauss])
        updateSliderInput(session, "W", value = GW[currGauss])
      })
      
      
      # Erneuert slider f?r den aktuellen Gauss (falls numGauss==0 --> lehrer output --> slider wird nicht angezeigt)
      observe({
        befehl$updateSliderCurrGauss
        if (numGauss>1){
          output$sliderCurrGauss <- renderUI({
            sliderInput('currGauss', h6('Gaussian No.'), min = 1, max = numGauss, value = currGauss, step= 1)
          })
          updateSliderInput(session, "currGauss", value = currGauss)
        } else {
          output$sliderCurrGauss <- renderUI({
            
          })
        }
      })
      #Erneuert Werte f?r M/S/W, wenn der Slider bewegt wird
      observe({
        #print("Refresh Values for M")
        if (!is.null(input$M)){
          GM[currGauss] <<- input$M
          iBefehl <<- iBefehl+1
          befehl$plot <- iBefehl
        }
      })
      observe({
        #print("Refresh Values for S")
        if (!is.null(input$S)){
          GS[currGauss] <<- input$S
          iBefehl <<- iBefehl+1
          befehl$plot <- iBefehl
        }
      })
      observe({
        #print("Refresh Values for W")
        if (!is.null(input$W)){
          GW[currGauss] <<- input$W
          iBefehl <<- iBefehl+1
          befehl$plot <- iBefehl
        }
      })
      
      observe({
        #print("Refresh numIterations")
        numIterations <<- input$numIterations
        if (numIterations<1) {numIterations <<- 1}
      })
      
      # Starte EMGauss()
      observe({ 
        GMSave <<- GM
        GSSave <<- GS
        GWSave <<- GW
        numGaussSave <<- numGauss 
        if (input$expMax>0){
          print("Expectation Maximation")
          Var <- EMGauss(data,GM,GS,GW,numIterations,fast=fast)
          GM <<- Var$Means
          GS <<- Var$SDs
          GW <<- Var$Weights
        }
        # L?sche Gauss mit Weight==0
        for (i in length(GW):1){
          if (GW[i]==0){
            if (i<numGauss){
              for (j in i:(numGauss-1)){
                GM[j] <<- GM[j+1]
                GS[j] <<- GS[j+1]
                GW[j] <<- GW[j+1]
              }
            }
            GM <<- GM[1:numGauss-1]
            GS <<- GS[1:numGauss-1]
            GW <<- GW[1:numGauss-1]   
            numGauss <<- numGauss-1 
            if (currGauss>i) currGauss <<- currGauss-1
            if (currGauss==i) currGauss <<- 1
          }
        }
        for (i in 1:numGauss){
          if (GM[i]<MLimit[1]) MLimit[1] <<- GM[i]*1.01^(sign(-GM[i]))   #Faktor 1.01^(sign(-GM[i]), weil nach dem runden von GM oder MLimit sonst Fehler auftreten k?nnen
          if (GM[i]>MLimit[2]) MLimit[2] <<- GM[i]*1.01^(sign(GM[i]))
          if (GS[i]<SLimit[1]) SLimit[1] <<- GS[i]*1.01^(sign(-GS[i]))
          if (GS[i]>SLimit[2]) SLimit[2] <<- GS[i]*1.01^(sign(GS[i]))
          if (GW[i]<WLimit[1]) WLimit[1] <<- GW[i]*1.01^(sign(-GW[i]))
          if (GW[i]>WLimit[2]) WLimit[2] <<- GW[i]*1.01^(sign(GW[i]))
        }
        iBefehl <<- iBefehl+1
        befehl$updateCurrGauss <- iBefehl
        befehl$updateSliderCurrGauss <- iBefehl
        befehl$drawSliderMSW <- iBefehl
        befehl$updateSlider <- iBefehl
        befehl$plot <- iBefehl
        befehl$updateMinMax <- iBefehl
      })
      
      # Lade Werte (wurden vor Expectation Maximation gespeichert)
      observe({
        if (input$restore>0){
          print("Restore previous Values")
          GM <<- GMSave
          GS <<- GSSave
          GW <<- GWSave
          numGauss <<- numGaussSave
          currGauss <<- 1
        }
        iBefehl <<- iBefehl+1
        befehl$updateSliderCurrGauss <- iBefehl
        befehl$updateSlider <- iBefehl
        befehl$plot <- iBefehl
      })
      
      #Lade Werte f?r best RMS
      observe({
        if (input$RestoreBestRMS>0){
          print("Restore Values of Best RMS")
          GM <<- GMBestRMS
          GS <<- GSBestRMS
          GW <<- GWBestRMS
          numGauss <<- numGaussBestRMS
          currGauss <<- 1
        }
        iBefehl <<- iBefehl+1
        befehl$updateSliderCurrGauss <- iBefehl
        befehl$updateSlider <- iBefehl
        befehl$plot <- iBefehl
      })
      
      # Plotten der Grafik (nur bei Befehl (befehl$plot))
      output$PDE <- renderPlot({
        befehl$plot
        #print("PDE estimation using ParetoDensityEstimation... ");
        #Plotte Pareto Density
        plot(Kernels, ParetoDensity,xlim=xlimit,ylim=ylimit, col="black", axes = FALSE, xlab = "Data", ylab = "Pareto Density Estimation", type="l", lwd=3,xaxs='i',yaxs='i')
        axis(1,xlim=xlimit,col="black",las=1) #x-Achse
        axis(2,ylim=ylimit,col="black",las=1) #y-Achse
        par(xaxs='i')
        par(yaxs='i')
        u <- par("usr")
        arrows(u[1], u[3], 1.05*u[2],u[3], code = 2, xpd = TRUE,lwd=2)
        arrows(u[1], u[3], u[1], 1.1*u[4], code = 2, xpd = TRUE,lwd=2)
        #box() #Kasten um Graphen
        par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
        # plotte gaussians (in Schwarz)
        FSum=0
        for (i in 1:numGauss){
          Fi=dnorm(Kernels,GM[i],GS[i])*GW[i]
          FSum=FSum+Fi
          if (input$showComponents){
            if (i==currGauss){
              points(Kernels, Fi,xlim=xlimit,ylim=ylimit, col="green", type="l", lwd=3)
            } else {
              points(Kernels, Fi,xlim=xlimit,ylim=ylimit, col="blue", type="l", lwd=3)
            }
            par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
          }
        }
        # Plotte Summe ?ber alle gaussians (in Rot)
        points(Kernels, FSum,xlim=xlimit,ylim=ylimit, col="red", type="l", lwd=3)
        par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
        #Plot Ende
        #Berechne RMS
        RMS <<- sqrt(sum((FSum-ParetoDensity)^2))/RMS0
        #output$RMS <- renderText({
         # paste("RMS = ",signif(RMS, digits=nSignif))
        #})
        str=paste('AdaptGauss(): GMM with',"RMS = ",signif(RMS, digits=nSignif))
        
        
        if (RMS<BestRMS){ # Speichere Werte, falls RMS optimiert wurde
          GMBestRMS <<- GM
          GSBestRMS <<- GS
          GWBestRMS <<- GW
          BestRMS <<- RMS
          numGaussBestRMS <<- numGauss
          currGaussBestRMS <<- currGauss
        }
        # Berechne (immer) und Plotte (nur wenn angeklickt) Bayes Boundaries
        if (numGauss>1 && sum(GW>0)>1){  # sum(GW>0)>1: Sicherstellen, dass mindestens 2 Gauss ein weight von ?ber 0 haben, sonst ist Berechnung von BayesBoundaries nicht m?glich
          BayesBoundaries <- BayesDecisionBoundaries(GM[GW>0],GS[GW>0],GW[GW>0])
          #BB <<- BayesBoundaries$DecisionBoundaries; 
          BB <<- BayesBoundaries
          if (input$showBayes){
            for (i in 1:length(BB)){
              abline(v=BB[i],col="Magenta")
              #plot(c(BB[i],BB[i]), ylimit,xlim=xlimit,ylim=ylimit, col="Magenta", axes = FALSE, xlab = " ", , ylab = " ", type="l", lwd=3)
              par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
            }
            BBText=paste(round(BB[1],digits=5-round(log(xlimit[2]-xlimit[1],10))))   #BB wird angepasst an den betrachteten X-Bereich (xlimit[2]-xlimit[1]) gerundet.
            if (length(BB)>1){
              for (i in 2:length(BB)){
                BBText=paste(BBText,", ",round(BB[i],digits=5-round(log(xlimit[2]-xlimit[1],10))))  #BB wird wieder angepasst an den betrachteten X-Bereich (xlimit[2]-xlimit[1]) gerundet.
              }
            }
            title(paste(str,", BayesBoundaries = ",BBText))
#             output$Bayes <- renderText({
#               paste("BayesBoundaries = ",BBText)
#             })
          } else {
            #output$Bayes <- renderText({  })
            title(str)
          }
        } 
        else{
          BB <<- NaN
          title(str)
          #output$Bayes <- renderText({  })
        }
      })
      
      # Fuege Gauss hinzu (bei Knopfdruck)
      observe({ 
        if (input$AddGaussButton > 0){
          print("Add Gauss")
          numGauss <<- numGauss+1
          GM[numGauss] <<- MDefault
          GS[numGauss] <<- SDefault
          GW[numGauss] <<- WDefault
          currGauss <<- numGauss
        }
        iBefehl <<- iBefehl+1
        befehl$updateSliderCurrGauss <- iBefehl
        befehl$plot <- iBefehl
      })
      
      # L?sche aktuellen Gauss: alle Gauss mit gr??erem Index rutschen eins nach link; neuer aktueller Gauss wird der, welcher vorher rechts vom gel?schten gauss stand. currGauss bleibt also gleich. Ausnahme: der Gaus mit dem gr??ten Index wird gel?scht, dann currGauss -> currGauss-1
      observe({
        if (input$DeleteGaussButton > 0 && numGauss>1){
          print("DeleteGauss")
          if (currGauss<numGauss){
            for (i in currGauss:(numGauss-1)){
              GM[i] <<- GM[i+1]
              GS[i] <<- GS[i+1]
              GW[i] <<- GW[i+1]
            }
          }
          if (currGauss==numGauss){
            currGauss <<- numGauss-1
          }
          GM <<- GM[1:numGauss-1]
          GS <<- GS[1:numGauss-1]
          GW <<- GW[1:numGauss-1]     
          numGauss <<- numGauss-1
        }
        iBefehl <<- iBefehl+1
        befehl$updateCurrGauss <- iBefehl
        befehl$updateSlider <- iBefehl
        befehl$updateSliderCurrGauss <- iBefehl
        befehl$plot <- iBefehl
      })
      
      
      #Closing APP oder plot figure
      observe({
        input$CloseButton
        input$PlotFig
        if (input$CloseButton+input$PlotFig > 0){
          print("Plot Figure")
          # Plotte Ergebnis
          plot(Kernels, ParetoDensity,xlim=xlimit,ylim=ylimit, col="black", axes = FALSE, xlab = "Data", ylab = "Pareto Density Estimation", main='GMM with AdaptGauss()',type="l", lwd=3)
          axis(1,xlim=xlimit,col="black",las=1) #x-Achse
          axis(2,ylim=ylimit,col="black",las=1) #y-Achse
          u <- par("usr")
          arrows(u[1], u[3], 1.05*u[2],u[3], code = 2, xpd = TRUE,lwd=2)
          arrows(u[1], u[3], u[1], 1.1*u[4], code = 2, xpd = TRUE,lwd=2)
          #box() #Kasten um Graphen
          par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
          FSum=0
          for (i in 1:numGauss){
            Fi=dnorm(Kernels,GM[i],GS[i])*GW[i]
            FSum=FSum+Fi
            if (isolate(input$showComponents)){
              points(Kernels, Fi,xlim=xlimit,ylim=ylimit, col="blue", type="l", lwd=3)
              par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
            }
          }
          
          if (input$showBayes && numGauss>1){
            for (i in 1:length(BB)){
              abline(v=BB[i], col="Magenta")
              #plot(c(BB[i],BB[i]), ylimit,xlim=xlimit,ylim=ylimit, col="Magenta", axes = FALSE, xlab = " ", , ylab = " ", type="l", lwd=3)
              par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
            }
          } 
          
          #Fi=dnorm(Kernels,GM[currGauss],GS[currGauss])*GW[currGauss]
          #FSum=FSum+Fi
          #plot(Kernels, Fi,xlim=xlimit,ylim=ylimit, col="green", axes = FALSE, xlab = " ", , ylab = " ", type="l", lwd=3)
          #par(new = TRUE) # der befehl das der n?chste den alten nicht ?bermalt
          plot(Kernels, FSum,xlim=xlimit,ylim=ylimit, col="red", axes = FALSE, xlab = " ", , ylab = " ", type="l", lwd=3)
          # Plot Ende
        }
        if (input$CloseButton > 0){
          print("close App")
          output <- list(Means=GM,SDs=GS,Weights=GW,ParetoRadius=ParetoRadius,RMS=RMS,BayesBoundaries=BB)
          stopApp(output)
        }
      })# end observe CloseButton and PlotFig
      
      session$onSessionEnded(function() {
        print("close App")
        output <- list(Means=GM,SDs=GS,Weights=GW,ParetoRadius=ParetoRadius,RMS=RMS,BayesBoundaries=BB)
        stopApp(output)
        # write out everything into files
      })
      
    }# end function server
    
    ))
  
  #outputApp=runApp(paste0(dbtDirectory(),'/dbt/dbt.EMforGauss/R/AdaptGauss2'));

  return(outputApp)

}


## Benutzte Funktionen: Subversion\PUB\dbt\...
# dbt.EMforGauss\EMGauss.R
# dbt.ClusteringAlgorithms\KmeansCluster.R
# dbt.BayesDecision\BayesDecisionBoundaries.R