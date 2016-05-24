require(shiny)
library(AFM)
library(rgl)
library(rglwidget)
library(tools)
library(data.table)
library(xtable)
library(ggplot2)
library(plyr)
library(scales)
library(fractaldim)
library(stringr)
library(grid)
library(gridExtra)
library(reshape2)
library(gstat)


data("AFMImageOfAluminiumInterface")
data("AFMImageOfNormallyDistributedHeights")
data("AFMImageOfOnePeak")
data("AFMImageOfRegularPeaks")

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 900MB.
options(shiny.maxRequestSize = 900*1024^2)


#HEADLESS<-TRUE
HEADLESS<-FALSE

testHEADLESS<-function() {
  
  result = tryCatch({
    open3d()
  }, warning = function(w) {
    print("warning rgl.open()")
    HEADLESS<-TRUE
  }, error = function(e) {
    print("error rgl.open()")
    HEADLESS<-TRUE
  }, finally = {
    if (!HEADLESS) rgl.close()
    print("finally rgl.open()")
  })
  
}

testHEADLESS()
#HEADLESS<-TRUE
print(paste("HEADLESS is", HEADLESS))
print(paste("interactive() is", interactive()))


#Do not open rgl windows with headless shiny server
if (HEADLESS) {
  options(rgl.useNULL = TRUE)
  
  devicesOpenedLength<-length(as.integer(rgl.dev.list()))
  while(devicesOpenedLength >0) {
    rgl.close()
    devicesOpenedLength<-devicesOpenedLength-1
  }
  #rgl.open()
}else{
  options(rgl.useNULL = FALSE)
}


afm_data_sets<-c("AFMImageOfAluminiumInterface","AFMImageOfNormallyDistributedHeights","AFMImageOfOnePeak","AFMImageOfRegularPeaks")

BrowserCanvasPixelLimit<-128


shinyServer(function(input, output, session) {
  
  shinyjs::disable("export3DModel3DButton")
  
  if (HEADLESS) {
    open3d()
    dev <- rgl.cur()
    save <- options(rgl.inShiny = TRUE)
    on.exit(options(save))
    
    session$onSessionEnded(function() {
      rgl.set(dev)
      rgl.close()
    })
    
    shinyjs::disable("snapshot3DButton")
    shinyjs::disable("displayIn3DFileButton")
    
    
    session$sendCustomMessage("sceneChange",
                              sceneChange("thewidget", skipRedraw = FALSE))
    session$onFlushed(function()
      session$sendCustomMessage("sceneChange",
                                sceneChange("thewidget", skipRedraw = FALSE)))
    
  }
  
  v <- reactiveValues(
    AFMImageAnalyser = NULL)
  
  clearData<-function() {
    disableButtons()
    v$AFMImageAnalyser<-NULL
  }
  
  disableButtons<-function() {
    shinyjs::disable("saveRdataFileButton")
    shinyjs::disable("RoughnessByLengthScaleButton")
    shinyjs::disable("downloadPSDPSDButton")
    shinyjs::disable("downloadRoughnessVsLengthscalePSDButton")
    shinyjs::disable("checkNormalityIsotropyCheckButton")
    shinyjs::disable("fitVariogramVarianceModelsButton")
    shinyjs::disable("calculateFractalDimensionsButton")
    shinyjs::disable("calculateNetworksNetworksButton")
    
    shinyjs::disable("displayIn3D3DButton")
    shinyjs::disable("snapshot3DButton")
    shinyjs::disable("calculate3DModel3DButton")
    shinyjs::disable("export3DModel3DButton")
    shinyjs::disable("generateCheckReport")
    shinyjs::disable("generateReport")
  }
  
  enableButtons<-function() {
    shinyjs::enable("saveRdataFileButton")
    shinyjs::enable("RoughnessByLengthScaleButton")
    shinyjs::enable("downloadPSDPSDButton")
    shinyjs::enable("downloadRoughnessVsLengthscalePSDButton")
    shinyjs::enable("checkNormalityIsotropyCheckButton")
    shinyjs::enable("fitVariogramVarianceModelsButton")
    shinyjs::enable("calculateFractalDimensionsButton")
    shinyjs::enable("calculateNetworksNetworksButton")
    
    shinyjs::enable("displayIn3D3DButton")
    if (!HEADLESS) shinyjs::enable("snapshot3DButton")
    shinyjs::enable("calculate3DModel3DButton")
    #shinyjs::enable("export3DModel3DButton")
    #shinyjs::enable("generateCheckReport")
    #shinyjs::enable("generateReport")
  }
  
  displayImageName<-function() {
    if (is.null(v$AFMImageAnalyser)) {
      print("v$AFMImageAnalyser null")
      return(NULL)
    }
    return(renderUI(HTML(c(paste0("<h4>Image</h4>",basename(v$AFMImageAnalyser@AFMImage@fullfilename), sep="")))))
  }
  
  #
  # Import/Export Data
  # 
  output$choose_inputtype <- renderUI({
    radioButtons("inputtype", "From", as.list(c("file", "dataset")), inline=TRUE)
  })
  
  # Check boxes
  output$choose_type <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(input$inputtype))
      return()
    clearData()
    if (input$inputtype=="file") {
      fileInput('file1', 'Load', accept=c('.txt','.rdata','.rda'))
    }else{
      radioButtons("choose_dataset", "Image", as.list(afm_data_sets))
    }
  })
  
  observeEvent(input$file1, {
    inFile <- input$file1
    clearData()
    enableButtons()
    if (file_ext(inFile$name)=="txt") {
      AFMImage<-importFromNanoscope(inFile$datapath)
      v$AFMImageAnalyser<-AFMImageAnalyser(copy(AFMImage))
      v$AFMImageAnalyser@AFMImage@fullfilename<-inFile$name
      #v$localfullfilename <-inFile$name
    }else{
      #print(inFile)
      x<-load(file= inFile$datapath)
      isolate({
        v$AFMImageAnalyser<-get(x)
        print(v$AFMImageAnalyser@fullfilename)
        #v$localfullfilename<-v$AFMImageAnalyser@fullfilename
        #v$AFMImageAnalyser@AFMImage<-v$AFMImageAnalyser@AFMImage
      })
      rm(x)
    }
    print("file1 button pushed")
  })
  
  observeEvent(input$choose_dataset, {
    inFile <- input$choose_dataset
    #print(inFile)
    if(is.null(inFile)) return()
    clearData()
    enableButtons()
    AFMImage<-get(inFile)
    v$AFMImageAnalyser<-AFMImageAnalyser(copy(AFMImage))
    #v$localfullfilename <-v$AFMImageAnalyser@AFMImage@fullfilename
    #print(v$localfullfilename)
    print("choose_dataset button pushed")
  })
  
  
  output$displayIn3DFileButton <- renderUI({
    if (HEADLESS == FALSE) {
      actionButton('displayIn3DFileButton', label = 'Display 3D model')
    }
  })
  
  observeEvent(input$displayIn3DFileButton, {
    if (is.null(input$displayIn3DFileButton)) return(NULL)
    #print(input$displayIn3DFileButton)
    if (input$displayIn3DFileButton==c(0)) return(NULL)
    #print(input$displayIn3DFileButton)
    if(is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@AFMImage))  return(NULL)
    displayIn3D(v$AFMImageAnalyser@AFMImage, 1024)
    print("displayIn3DFileButton button pushed")
  })
  
  output$saveRdataFileButton <- downloadHandler(
    filename = function() { paste(basename(v$AFMImageAnalyser@AFMImage@fullfilename), '-AFMImageAnalyser','.Rda', sep='') },
    content = function(file) {
      if(!is.null(v$AFMImageAnalyser)) {
        print("Exporting calculation")
        AFMImageAnalyser=copy(v$AFMImageAnalyser)
        save(AFMImageAnalyser, file= file)      
        print("done")
      }
    }
  )
  
  #
  # PSD tab observer
  #
  output$downloadPSDPSDButton <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@psdAnalysis))
      return(NULL)
    downloadButton('exportPSD',label='Export PSD')
  })
  
  output$exportPSD <- downloadHandler(
    filename = function() { paste(basename(v$AFMImageAnalyser@AFMImage@fullfilename), '.csv', sep='') },
    content = function(file) {
      write.csv(v$AFMImageAnalyser@psdAnalysis@psd1d, file, row.names = FALSE)
    }
  )
  
  
  output$downloadRoughnessVsLengthscalePSDButton <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@psdAnalysis))
      return(NULL)
    downloadButton('exportRoughnessVsLengthscale',label='Export roughness vs. lengthscale')
  })
  output$exportRoughnessVsLengthscale <- downloadHandler(
    filename = function() { paste(basename(v$AFMImageAnalyser@AFMImage@fullfilename), '.csv', sep='') },
    content = function(file) {
      write.csv(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale, file, row.names = FALSE)
    }
  )
  
  
  observeEvent(input$RoughnessByLengthScaleButton, {
    input$RoughnessByLengthScaleButton
    print("RoughnessByLengthScaleButton button pushed")
    if (is.null(input$RoughnessByLengthScaleButton)) {
      print("input$RoughnessByLengthScaleButton==NULL")
      return(NULL)
    }
    print("input$RoughnessByLengthScaleButton!=NULL")
    if(input$RoughnessByLengthScaleButton == c(0))
    {
      print("input$RoughnessByLengthScaleButton==0")
      return()
    }else{
      isolate({
        input$RoughnessByLengthScaleButton
        
        # Create a Progress object
        progressPSD <- shiny::Progress$new()
        #progressPSD$set(message = "Calculting", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progressPSD$close())
        
        print("calculation of PSD")
        print(paste("with", 2^input$breaksSliderPSD, "breaks"))
        
        #createAFMImageAnalyser()
        psdAnalysis<-AFMImagePSDAnalysis()
        # Create a closure to update progress
        psdAnalysis@updateProgress<- function(value = NULL, detail = NULL, message = NULL) {
          if (exists("progressPSD")){
            if (!is.null(message)) {
              progressPSD$set(message = message, value = 0)
            }else{
              progressPSD$set(value = value, detail = detail)
            }
          }
        }
        psdAnalysis@psd1d_breaks<-2^input$breaksSliderPSD
        psdAnalysis@psd2d_truncHighLengthScale<-TRUE
        psdAnalysis<-performAllPSDCalculation(AFMImagePSDAnalysis= psdAnalysis, AFMImage= v$AFMImageAnalyser@AFMImage)
        print("done psdAnalysis")
        
        v$AFMImageAnalyser@psdAnalysis<-psdAnalysis
        print("done v$AFMImageAnalyser@psdAnalysis<-psdAnalysis")
      })
    }
    
    
  })
  
  #
  # Variance checks tab
  #
  output$imageNameCheck<-renderUI({
    imageName<-displayImageName()
    
    if (is.null(imageName)) {
      output$imageNameCheck<-renderUI(HTML(c("<h4>please select image first</h4>")))
      return(NULL)
    }
    
    output$imageNameCheck<-imageName
    print(imageName)
  })
  
  output$normalityVarianceCheckUI<-renderUI({
    input$checkNormalityIsotropyCheckButton
    normalityVarianceCheckUI()
  })
  normalityVarianceCheckUI<-reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@variogramAnalysis)) {
      return(list())
    }
    return(h3("Normality check"))
  })
  
  output$isotropyVarianceCheckUI<-renderUI({
    input$checkNormalityIsotropyCheckButton
    isotropyVarianceCheckUI()
  })
  isotropyVarianceCheckUI<-reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@variogramAnalysis)) {
      return(list())
    }
    return(h3("Isotropy checks"))
  })
  
  
  
  output$normalityIsotropyVarianceCheckImage<- renderImage({
    normalityIsotropyVarianceCheckImage()
  }, deleteFile = TRUE)
  normalityIsotropyVarianceCheckImage <- reactive({
    input$checkNormalityIsotropyCheckButton
    print("checkNormalityIsotropyCheckButton button pushed")
    
    if (is.null(input$checkNormalityIsotropyCheckButton)) {
      print("input$checkNormalityIsotropyCheckButton==NULL")
      return(NULL)
    }
    print("input$checkNormalityIsotropyCheckButton!=NULL")
    
    if (is.null(v$AFMImageAnalyser)||
        is.null(v$AFMImageAnalyser@variogramAnalysis)||
        length(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms)==0
        #||        input$checkNormalityIsotropyCheckButton == c(0)
    )  {
      return(list(src = tempfile()))
    }else{
      outfile1 <- tempfile(fileext='.png')
      png(outfile1, width=640, height=400)
      checkNormality(AFMImage= v$AFMImageAnalyser@AFMImage, v$AFMImageAnalyser@variogramAnalysis)
      dev.off()
      
      return(list(src = outfile1,
                  contentType = 'image/png',
                  width = 640,
                  height = 400,
                  alt = "fd2d_squareincr"))
      
    }
    
    return(list(src = tempfile()))
  })
  
  output$directionalVariogramsVarianceCheckImage<- renderPlot({
    directionalVariogramsVarianceCheckPlot()
  })
  directionalVariogramsVarianceCheckPlot <- reactive({
    input$checkNormalityIsotropyCheckButton
    print("checkNormalityIsotropyCheckButton button pushed")
    
    if (is.null(input$checkNormalityIsotropyCheckButton)) {
      print("input$checkNormalityIsotropyCheckButton==NULL")
      return(NULL)
    }
    print("input$checkNormalityIsotropyCheckButton!=NULL")
    
    if (is.null(v$AFMImageAnalyser)||
        is.null(v$AFMImageAnalyser@variogramAnalysis)||
        is.null(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms)||
        length(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms)==0) {
      #|| input$checkNormalityIsotropyCheckButton == c(0)
      return(NULL)
    }else{
      
      p2 <- ggplot(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms, 
                   aes(x=dist, y=gamma, color= as.factor(dir.hor), shape=as.factor(dir.hor)))
      p2 <- p2 + expand_limits(y = 0)
      p2 <- p2 + geom_point()
      p2 <- p2 + geom_line()
      p2 <- p2 + ylab("semivariance (nm^2)")
      p2 <- p2 + xlab("distance (nm)")
      p2 <- p2 + ggtitle("Directional variograms")
      return(p2)
    }
    return(NULL)
  })
  
  observeEvent(input$checkNormalityIsotropyCheckButton, {
    input$checkNormalityIsotropyCheckButton
    print("checkNormalityIsotropyCheckButton button pushed")
    
    if (is.null(input$checkNormalityIsotropyCheckButton)) {
      print("input$checkNormalityIsotropyCheckButton==NULL")
      return(NULL)
    }
    print("input$checkNormalityIsotropyCheckButton!=NULL")
    
    if(input$checkNormalityIsotropyCheckButton == c(0)) {
      print("input$checkNormalityIsotropyCheckButton==c(0)")
    }else{
      isolate({
        # Create a Progress object
        progressVariogramAnalysis <- shiny::Progress$new()
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progressVariogramAnalysis$close())
        # Create a closure to update progress.
        print("Calculation of directional variograms")
        #createAFMImageAnalyser()
        sampleFitPercentage<-3.43/100
        variogramAnalysis<-AFMImageVariogramAnalysis(sampleFitPercentage= sampleFitPercentage)
        variogramAnalysis@updateProgress<-function(value = NULL, detail = NULL, message = NULL) {
          if (!is.null(message)) {
            progressVariogramAnalysis$set(message = message, value = 0)
          }else{
            progressVariogramAnalysis$set(value = value, detail = detail)
          }
        }
        variogramAnalysis@updateProgress(message="Calculating directional variograms", value=0)
        variogramAnalysis@updateProgress(value= 0, detail = "0/1")
        variogramAnalysis@directionalVariograms<- calculateDirectionalVariograms(AFMImage= sampleAFMImage(v$AFMImageAnalyser@AFMImage,input$sampleIsotropyVarianceCheckSlider),AFMImageVariogramAnalysis= variogramAnalysis)
        
        # if models were already calculated, do not erase them
        if(!is.null(v$AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram)&&
           length(v$AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram)!=0) {
          v$AFMImageAnalyser@variogramAnalysis@directionalVariograms<-variogramAnalysis@directionalVariograms
        }else{
          v$AFMImageAnalyser@variogramAnalysis<-variogramAnalysis
        }
        
        
        
        
        print("done")
      })
    }
    return(NULL)  
  })
  
  #
  # Variance model tab
  #
  output$imageNameVarianceModels<-renderUI({
    imageName<-displayImageName()
    if (is.null(imageName)) {
      output$imageNameVarianceModels<-renderUI(HTML(c("<h4>please select image first</h4>")))
      return(NULL)
    }
    output$imageNameVarianceModels<-imageName
    print(imageName)
  })
  
  output$bestmodeltableVarianceModelsUI<-renderUI({
    input$fitVariogramVarianceModelsButton
    bestmodeltableVarianceModelsUI()
  })
  bestmodeltableVarianceModelsUI<-reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@variogramAnalysis@variogramModels)) {
      return(list())
    }
    return(h3("Variogram models"))
  })
  
  
  output$bestmodeltableVarianceModelsPlot <- renderTable({ 
    input$fitVariogramVarianceModelsButton
    bestmodeltableVarianceModelsPlot()
  }, include.rownames=FALSE, include.colnames=TRUE)
  bestmodeltableVarianceModelsPlot <- reactive({
    input$fitVariogramVarianceModelsButton
    print("fitVariogramVarianceModelsButton button pushed")
    
    if (is.null(input$fitVariogramVarianceModelsButton)) {
      print("input$fitVariogramVarianceModelsButton==NULL")
      return(NULL)
    }
    print("input$fitVariogramVarianceModelsButton!=NULL")
    #     if(input$fitVariogramVarianceModelsButton == c(0)) {
    #       print("input$fitVariogramVarianceModelsButton==c(0)")
    #       return(NULL)
    #     }   
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@variogramAnalysis@variogramModels)) {
      return(NULL)
    }
    
    # get variogram model evaluation
    if (length(v$AFMImageAnalyser@variogramAnalysis@variogramModels)!=0) {
      
      mergedDT<-getDTModelEvaluation(v$AFMImageAnalyser@variogramAnalysis)
      #print(mergedDT)
      
      sillrangeDT<-getDTModelSillRange(v$AFMImageAnalyser@variogramAnalysis)
      setkey(sillrangeDT, "model")
      name<-press<-NULL
      sampleDT <- mergedDT[name==basename(v$AFMImageAnalyser@AFMImage@fullfilename)]
      setkey(sampleDT, "model")
      #sampleDT <- sampleDT[cor>0.98]
      sampleDT<-merge(sampleDT, sillrangeDT, by="model")
      sampleDT<-sampleDT[,name:=NULL]
      sampleDT <- unique(sampleDT)
      sampleDT <- sampleDT[order(-rank(cor), rank(press))]
      #print(sampleDT)
      
      tot<-sum(v$AFMImageAnalyser@AFMImage@data$h)
      
      summarySampleDT<-copy(sampleDT)
      summarySampleDT$press<-round(sampleDT$press)
      summarySampleDT$sill<-round(sampleDT$sill)
      summarySampleDT$range<-round(sampleDT$range)
      
      print("plotting variogram table...")
      
      fractalDF<-data.frame(summarySampleDT)
      #print(fractalDF)
      return({xtable(fractalDF)})
    }else{
      return(NULL)
    }
  })
  
  getSpplotColors<-function(colLimit) {
    blues9[3:colLimit]  
  }
  
  output$allmodelsModelImage<-renderImage({
    allmodelsModelImage()
  })
  allmodelsModelImage<- reactive({
    input$fitVariogramVarianceModelsButton
    print("fitVariogramVarianceModelsButton button pushed")
    
    if (is.null(input$fitVariogramVarianceModelsButton)) {
      print("input$fitVariogramVarianceModelsButton==NULL")
      outfile1 <- tempfile(fileext='.png')
      png(outfile1, width=800, height=300)
      dev.off()
      return(list(src = outfile1,
                  contentType = 'image/png',
                  width = 800,
                  height = 300,
                  alt = "Models not calculated"))
    }
    print("input$fitVariogramVarianceModelsButton!=NULL")
    
    #     if(input$fitVariogramVarianceModelsButton == c(0)) {
    #       print("input$fitVariogramVarianceModelsButton==c(0)")
    #       outfile1 <- tempfile(fileext='.png')
    #       png(outfile1, width=800, height=300)
    #       dev.off()
    #       return(list(src = outfile1,
    #                   contentType = 'image/png',
    #                   width = 800,
    #                   height = 300,
    #                   alt = "Models not calculated"))
    #     }
    if (is.null(v$AFMImageAnalyser)||
        is.null(v$AFMImageAnalyser@variogramAnalysis)||
        is.null(v$AFMImageAnalyser@variogramAnalysis@variogramModels)) {
      outfile1 <- tempfile(fileext='.png')
      png(outfile1, width=800, height=300)
      dev.off()
      return(list(src = outfile1,
                  contentType = 'image/png',
                  width = 800,
                  height = 300,
                  alt = "Models not calculated"))
    }
    
    # get variogram model evaluation
    if (length(v$AFMImageAnalyser@variogramAnalysis@variogramModels)!=0) {
      mergedDT<-getDTModelEvaluation(v$AFMImageAnalyser@variogramAnalysis)
      #print(mergedDT)
      
      sillrangeDT<-getDTModelSillRange(v$AFMImageAnalyser@variogramAnalysis)
      setkey(sillrangeDT, "model")
      name<-press<-NULL
      sampleDT <- mergedDT[name==basename(v$AFMImageAnalyser@AFMImage@fullfilename)]
      setkey(sampleDT, "model")
      #sampleDT <- sampleDT[cor>0.98]
      sampleDT<-merge(sampleDT, sillrangeDT, by="model")
      sampleDT<-sampleDT[,name:=NULL]
      sampleDT <- unique(sampleDT)
      sampleDT <- sampleDT[order(-rank(cor), rank(press))]
      
      #####################
      # new page for experimental variogram and models
      numberOfModelsPerPage=length(v$AFMImageAnalyser@variogramAnalysis@variogramModels)
      
      allVarioModels<-str_sub(sampleDT$model,-3)
      
      outfile1 <- tempfile(fileext='.png')
      #print(outfile1)
      png(outfile1, width=800, height=numberOfModelsPerPage*300)
      
      #grid.newpage()
      printVariogramModelEvaluations(AFMImageAnalyser= v$AFMImageAnalyser, sampleDT= sampleDT, numberOfModelsPerPage= numberOfModelsPerPage)
      dev.off()
      
      return(list(src = outfile1,
                  contentType = 'image/png',
                  width = 800,
                  height = 1600,
                  alt = "all evaluated models"))
    }
    
  })
  
  observeEvent(input$fitVariogramVarianceModelsButton, {
    input$fitVariogramVarianceModelsButton
    print("fitVariogramVarianceModelsButton button pushed")
    
    if (is.null(input$fitVariogramVarianceModelsButton)) {
      print("input$fitVariogramVarianceModelsButton==NULL")
      return(NULL)
    }
    print("input$fitVariogramVarianceModelsButton!=NULL")
    
    if(input$fitVariogramVarianceModelsButton == c(0)) {
      print("input$fitVariogramVarianceModelsButton==c(0)")
    }else{
      isolate({
        # Create a Progress object
        progressVariogramAnalysis <- shiny::Progress$new()
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progressVariogramAnalysis$close())
        
        
        #createAFMImageAnalyser()
        sampleFitPercentage<-input$sampleFitVarianceModelsSlider/100
        
        
        variogramAnalysis<-AFMImageVariogramAnalysis(sampleFitPercentage= sampleFitPercentage)
        variogramAnalysis@updateProgress<-function(value = NULL, detail = NULL, message = NULL)  {
          if (is.null(value)&& is.null(detail)&& is.null(message)) return(TRUE)
          if (!is.null(message)) {
            progressVariogramAnalysis$set(message = message, value = 0)
          }else{
            progressVariogramAnalysis$set(value = value, detail = detail)
          }
        }
        
        variogramAnalysis@updateProgress(message="Calculating omnidirectional variogram", value=0)
        variogramAnalysis@omnidirectionalVariogram<- calculateOmnidirectionalVariogram(AFMImage=v$AFMImageAnalyser@AFMImage,AFMImageVariogramAnalysis= variogramAnalysis)
        variogramAnalysis@updateProgress(message="Fitting and evaluating",value=0)
        
        # normality and isotropy already calculated
        if(!is.null(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms)&&
           length(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms)!=0) {
          variogramAnalysis@directionalVariograms<-v$AFMImageAnalyser@variogramAnalysis@directionalVariograms
        }
        v$AFMImageAnalyser@variogramAnalysis<-evaluateVariogramModels(variogramAnalysis, v$AFMImageAnalyser@AFMImage)
        print("done")
      })
    }
    return(NULL)  
  })
  #
  # File Tab display
  #
  output$basicInfoFileTable <- renderTable({ imageInformations()}, include.rownames=FALSE, include.colnames=FALSE)
  output$roughnessesFileTable <- renderTable({ roughnessesImagesTable()}, include.rownames=FALSE, include.colnames=FALSE)
  
  output$imageInformationsUI<-renderUI({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@AFMImage)) return(NULL)
    h3("Image informations")
  })
  output$roughnessUI<-renderUI({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@AFMImage)) return(NULL)
    h3("Roughnesses")
  })
  
  
  imageInformations <- function(){ 
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@AFMImage)) return(NULL)
    
    scansize<-v$AFMImageAnalyser@AFMImage@scansize
    samplesperline<-as.character(v$AFMImageAnalyser@AFMImage@samplesperline)
    lines<-as.character(v$AFMImageAnalyser@AFMImage@lines)
    
    charact<-data.frame(name=c("scansize","samplesperline","lines"),values=c(paste0(scansize,"nm"),paste0(samplesperline,"px"),paste0(lines,"px")))
    return({xtable(charact)})
  }
  
  roughnessesImagesTable <- function(){ 
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@AFMImage)) return(NULL)
    roughnesses<-getRoughnessParameters(v$AFMImageAnalyser@AFMImage)
    roughnessDF<-data.frame(name=c("total Rrms","Ra (mean roughness)"),values=c(paste0(round(roughnesses$totalRMSRoughness_TotalRrms, digits=4),"nm"),paste0(round(roughnesses$MeanRoughness_Ra, digits=4),"nm")))
    return({xtable(roughnessDF)})
  }
  
  #
  # PSD Tab display
  #
  output$imageNamePSD<-renderUI({
    imageName<-displayImageName()
    
    if (is.null(imageName)) {
      output$imageNamePSD<-renderUI(HTML(c("<h4>please select image first</h4>")))
      return(NULL)
    }
    
    output$imageNamePSD<-imageName
    #print(imageName)
  })
  
  output$plotPSDUI<-renderUI({
    myPlotPSDUI()
  })
  
  myPlotPSDUI<-reactive({
    if (is.null(v$AFMImageAnalyser)||
        is.null(v$AFMImageAnalyser@psdAnalysis)) {
      return(NULL)
    }
    h3("Power spectrum density")
  })
  
  
  output$plotPSDRvsLUI<-renderUI({
    myplotPSDRvsLUI()
  })
  
  myplotPSDRvsLUI<-reactive({
    if (is.null(v$AFMImageAnalyser)||
        is.null(v$AFMImageAnalyser@psdAnalysis)) {
      return(NULL)
    }
    h3("Roughness vs. lengthscale")
  })
  
  
  output$plotPSD <- renderPlot({
    my_PSD_plot()
  })
  
  my_PSD_plot <- reactive({
    if (is.null(input$RoughnessByLengthScaleButton)) {
      print("input$RoughnessByLengthScaleButton==NULL")
      return(NULL)
    }
    print("input$RoughnessByLengthScaleButton!=NULL")
    
    if (length(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)==0) {
      print("length(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)==0")
      return(NULL)
    }
    
    #     if(input$RoughnessByLengthScaleButton == c(0)) {
    #       print("my_PSD_plot input$RoughnessByLengthScaleButton==0")
    #     }else{
    isolate({
      if (is.null(v$AFMImageAnalyser)) {
        print("is.null(v$AFMImageAnalyser)")
      }else{ if (is.null(v$AFMImageAnalyser@psdAnalysis)) {
        print("is.null(v$AFMImageAnalyser@psdAnalysis)")
      }else{
        datap<-v$AFMImageAnalyser@psdAnalysis@psd1d
        p <- ggplot(data=datap)
        p <- p + geom_point(aes(freq, PSD, color=type),data=datap[datap$type %in% c("PSD-2D")])
        p <- p + geom_line(aes(freq, PSD, color=type),data=datap[datap$type %in% c("PSD-1D")],size=1.1)
        p <- p + scale_x_log10()
        p <- p + scale_y_log10()
        p <- p + ylab("PSD (nm^4)")
        p <- p + xlab("Frequency (nm^-1)")
        return(p)
      }}
    })
    # }
    return(NULL)
  })
  
  
  output$plotPSDRvsL <- renderPlot({
    my_RoughnessVsLengthscale_plot()
  })
  
  my_RoughnessVsLengthscale_plot <- reactive({
    if (is.null(input$RoughnessByLengthScaleButton)) {
      print("input$RoughnessByLengthScaleButton==NULL")
      return(NULL)
    }
    print("input$RoughnessByLengthScaleButton!=NULL")
    print(length(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale))
    if (length(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)==0) {
      print("length(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)==0")
      return(NULL)
    }
    
    #     if(input$RoughnessByLengthScaleButton == c(0))
    #     {
    #       print("my_RoughnessVsLengthscale_plot input$RoughnessByLengthScaleButton==0")
    #     }else{
    isolate({
      input$RoughnessByLengthScaleButton
      
      if (is.null(v$AFMImageAnalyser)) {
        print("is.null(v$AFMImageAnalyser)")
      }else{ if (is.null(v$AFMImageAnalyser@psdAnalysis)) {
        print("is.null(v$AFMImageAnalyser@psdAnalysis)")
      }else{
        print("Displaying Roughness vs. lengthscale")
        r<-roughness<-filename<-NULL
        #v$AFMImagePSDAnalysis@roughnessAgainstLengthscale$filename<-v$localfullfilename
        p1 <- ggplot(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale, aes(x=r, y=roughness, colour= "red"))
        p1 <- p1 + geom_point()
        p1 <- p1 + geom_line()
        p1 <- p1 + ylab("roughness (nm)")
        p1 <- p1 + xlab("lengthscale (nm)")
        p1 <- p1 + guides(color=FALSE)
        p1 <- p1 + scale_color_discrete(NULL)
        return(p1)
      }}
    })
    # }
    return(NULL)
  })
  
  
  #
  # 3D display
  #
  output$imageName3D<-renderUI({
    imageName<-displayImageName()
    if (is.null(imageName)) {
      output$imageName3D<-renderUI(HTML(c("<h4>please select image first</h4>")))
      return(NULL)
    }
    print(imageName)
    output$imageName3D<-imageName 
  })
  
  output$panel3DUI<-renderUI({
    input$calculateFractalDimensionsButton
    myPanel3DUI()
  })
  myPanel3DUI<-reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@AFMImage)) {
      return(list())
    }
    if (!HEADLESS) return(h3("Three dimensional display"))
    
    if ((v$AFMImageAnalyser@AFMImage@samplesperline>BrowserCanvasPixelLimit)||
        (v$AFMImageAnalyser@AFMImage@lines>BrowserCanvasPixelLimit))
      return(h3("Simplified three dimensional display"))
    else
      return(h3("Three dimensional display"))
  })
  
  observeEvent(input$displayIn3D3DButton, {
    input$displayIn3D3DButton
    if (is.null(input$displayIn3D3DButton)) return(NULL)
    if (input$displayIn3D3DButton==c(0)) return(NULL)
    if(is.null(v$AFMImageAnalyser@AFMImage))  return(NULL)
    print("displayIn3D3DButton button pushed observeEvent")
    
    copyH<-copy(v$AFMImageAnalyser@AFMImage@data$h)
    
    width<-600
    #AFMImage<-v$AFMImageAnalyser@AFMImage
    if (HEADLESS) AFMImage<-simplifyAFMImage(v$AFMImageAnalyser@AFMImage, BrowserCanvasPixelLimit,BrowserCanvasPixelLimit)
    else AFMImage<-v$AFMImageAnalyser@AFMImage
    
    # respect the proportion between horizontal / vertical distance and heigth
    newHeights <- input$height3Dslider*(AFMImage@data$h)*(AFMImage@samplesperline)/(AFMImage@scansize)
    minH<-min(newHeights)
    # TODO check validity of created image instead
    if(!is.na(minH)) {
      newH<-(newHeights-minH)
      y<-matrix(newH, nrow = AFMImage@lines, ncol = AFMImage@samplesperline)
      z <- seq(1,ncol(y),by=1) 
      x <- (1:nrow(y))
      
      ylim <- range(y)
      ylen <- ylim[2] - ylim[1] + 1
      print(ylen)
      colorlut <- heat.colors(ylen, alpha = 1) # height color lookup table
      col <- colorlut[ y-ylim[1]+1 ] # assign colors to heights
      
      print(HEADLESS)
      if ((!HEADLESS)||is.null(rgl.cur())||rgl.cur()==0) { 
        print("open rgl")
        rgl.open()
      }else{
        print("using dev")
        dev <- rgl.cur()
        rgl.set(dev)
        rgl.clear()
      }
      
      bboxylen=3
      if(ylim[2]<60) bboxylen=2
      
      
      if (!HEADLESS) {
        par3d(windowRect = 100 + c( 0, 0, width, width ) )
        rgl.bg(color = c("white"),  back = "lines")
        
        rgl.bbox(color = c("#333333", "black"), emission = "#333333", 
                 specular = "#111111", shininess = 0, alpha = 0.6, xlen=0, zlen=0, ylen=bboxylen )
        rgl.surface(x, z, y, color=col, back="lines")
        i<-130
        rgl.viewpoint(i,i/4,zoom=1.1)
      }else{
        
        if (!exists("mysurface")||is.null(mysurface)) {
          myb3d<-bg3d(color = c("white"),  back = "lines")
          axes3d(color = c("#333333", "black"),
                 labels=FALSE)  
          #           mybbox3d<-bbox3d(color = c("#333333", "black"), 
          #                            emission = "#000000", 
          #                  #specular = "#111111",
          #                  specular = "#000000",
          #                  front="line", back="line",
          #                  shininess = 0, draw_front=FALSE, alpha = 0.6, 
          #                  xlen=0, zlen=0, ylen=bboxylen )
          
          s3d<-surface3d(x, z, y, color=col, back="lines")
          print(paste("add", s3d))
          root <- currentSubscene3d()
          newSubscene3d("inherit", "inherit", "inherit", copyShapes = TRUE, parent = root)
          clipplanes3d(1, 0, 0, 0)
          mysurface <- scene3d()
          plot3d(mysurface)        
        }else{
          print(paste("delete", s3d))
          useSubscene3d(root)
          delFromSubscene3d(s3d)
          s3d<-surface3d(x, z, y, color=col, back="lines")
        }
        i<-130
        rgl.viewpoint(i,i/4,zoom=1.1)
      }
      
    }
    v$AFMImageAnalyser@AFMImage@data$h<-copyH
    print("end display")
    # })
  })
  
  #   output$thewidget <- renderRglwidget({
  #         if (is.null(input$displayIn3D3DButton)) return(rglwidget())
  #         if (input$displayIn3D3DButton==c(0)) return(rglwidget())
  #         if(is.null(v$AFMImageAnalyser@AFMImage))  return(rglwidget())
  #     print("displayIn3D3DButton button pushed renderRglwidget")
  #     rglwidget()
  #     print("end renderRglwidget")
  #   })
  
  #     output$thecontroller <-
  #       renderRglcontroller({
  #   #       if (input$displayIn3D3DButton!=c(0))
  #   #         # It makes more sense to use rglcontroller as below, but
  #   #         # this works...
  #            playwidget("thewidget", respondTo = "height3Dslider")
  #       })  
  #   
  
  output$thewidget <- renderRglwidget({
    if (HEADLESS) {
      if (is.null(input$displayIn3D3DButton)) return(NULL)
      if (input$displayIn3D3DButton==c(0)) return(NULL)
      if(is.null(v$AFMImageAnalyser@AFMImage))  return(NULL)
      # print("displayIn3D3DButton button pushed renderRglwidget")    
      rglwidget(height=600, reuse=FALSE)
      #rglwidget(height=600, controllers="thecontroller")
    }
  })
  
  
  output$snapshot3DButton <- downloadHandler(
    filename = function() { paste(basename(v$AFMImageAnalyser@AFMImage@fullfilename), '-multiply', input$height3Dslider,'-3D.png', sep='') },
    content = function(file) {
      if(!is.null(v$AFMImageAnalyser@AFMImage)) {
        print("Exporting 3D png")
        rgl.snapshot(filename = file)
        print("done")
      }
    }
  )
  
  
  observeEvent(input$calculate3DModel3DButton, {
    input$calculate3DModel3DButton
    print("calculate3DModel3DButton button pushed")
    if (is.null(input$calculate3DModel3DButton)) {
      print("input$calculate3DModel3DButton==NULL")
      return(NULL)
    }
    print("input$calculate3DModel3DButton!=NULL")
    if(input$calculate3DModel3DButton == c(0))
    {
      print("input$calculate3DModel3DButton==0")
      return()
    }else{
      isolate({
        input$calculate3DModel3DButton
        
        # Create a Progress object
        progressCalculate3DModel <- shiny::Progress$new()
        #progressPSD$set(message = "Calculting", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progressCalculate3DModel$close())
        
        print("calculation of 3D Model")
        
        #createAFMImageAnalyser()
        modelAnalysis<-new ("AFMImage3DModelAnalysis")
        # Create a closure to update progress
        modelAnalysis@updateProgress<- function(value = NULL, detail = NULL, message = NULL) {
          if (!is.null(message)) {
            progressCalculate3DModel$set(message = message, value = 0)
          }else{
            progressCalculate3DModel$set(value = value, detail = detail)
          }
          return(TRUE)
        }
        modelAnalysis@updateProgress(message="Calculating 3D faces", value=0)
        
        if(!is.null(v$AFMImageAnalyser@AFMImage)) {
          print("Exporting 3D printing model")
          copyH<-copy(v$AFMImageAnalyser@AFMImage@data$h)
          v$AFMImageAnalyser@AFMImage@data$h<-v$AFMImageAnalyser@AFMImage@data$h*input$height3Dslider
          modelAna<-calculate3DModel(AFMImage3DModelAnalysis= modelAnalysis, AFMImage= v$AFMImageAnalyser@AFMImage)
          v$AFMImageAnalyser@AFMImage@data$h<-copyH
        }
        
        print("calculate3DModel done")
        
        v$AFMImageAnalyser@threeDimensionAnalysis<-modelAna
        print("done v$AFMImageAnalyser@threeDimensionAnalysis<-modelAnalysis")
        shinyjs::enable("export3DModel3DButton")
      })
    }
    
    
  })
  
  output$export3DModel3DButton <- downloadHandler(
    filename = function() { paste(basename(v$AFMImageAnalyser@AFMImage@fullfilename), '-multiply', input$height3Dslider, '.stl', sep='') },
    content = function(file) {
      if(!is.null(v$AFMImageAnalyser@AFMImage)) {
        print("Exporting 3D printing model")
        
        copyH<-copy(v$AFMImageAnalyser@AFMImage@data$h)
        v$AFMImageAnalyser@AFMImage@data$h<-v$AFMImageAnalyser@AFMImage@data$h*input$height3Dslider
        
        exportToSTL(AFMImage3DModelAnalysis=v$AFMImageAnalyser@threeDimensionAnalysis, 
                    AFMImage=v$AFMImageAnalyser@AFMImage, 
                    stlfullfilename=file)
        v$AFMImageAnalyser@AFMImage@data$h<-copyH
        print("done")
      }
    }
  )
  
  #
  # Fractal tab observer
  #
  observeEvent(input$calculateFractalDimensionsButton, {
    
    input$calculateFractalDimensionsButton
    print("calculateFractalDimensionsButton button pushed")
    
    print("input$calculateFractalDimensionsButton!=NULL")
    if(input$calculateFractalDimensionsButton == c(0))
    {
      print("input$calculateFractalDimensionsButton==0")
      return()
    }else{
      #      isolate({
      input$calculateFractalDimensionsButton
      
      # Create a Progress object
      progressFractal <- shiny::Progress$new()
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progressFractal$close())
      
      print("calculation of Fractal dimensions and scales")
      #createAFMImageAnalyser()
      # fractal dimension analysis
      fdAnalysis<-AFMImageFractalDimensionsAnalysis()
      # Create a closure to update progress
      fdAnalysis@updateProgress<-function(value = NULL, detail = NULL, message = NULL) {
        if (!is.null(message)) {
          progressFractal$set(message = message, value = 0)
        }else{
          progressFractal$set(value = value, detail = detail)
        }
      }
      fdAnalysis@fractalDimensionMethods<-getFractalDimensions(v$AFMImageAnalyser@AFMImage, fdAnalysis)
      
      print("done v$AFMImageAnalyser@fdAnalysis<-fdAnalysis")
      #      })
      v$AFMImageAnalyser@fdAnalysis<-fdAnalysis
    }
  })
  
  #
  # Fractal calculation display
  # 
  output$imageNameFractal<-renderUI({
    imageName<-displayImageName()
    if (is.null(imageName)) {
      output$imageNameFractal<-renderUI(HTML(c("<h4>please select image first</h4>")))
      return(NULL)
    }
    output$imageNameFractal<-imageName
    # print(imageName)
  })
  
  output$fractalDimensionsFractalTable <- renderTable({ 
    input$calculateFractalDimensionsButton
    fractalDimensionsFractalTable()
  }, include.rownames=FALSE, include.colnames=TRUE)
  
  fractalDimensionsFractalTable <- reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@fdAnalysis)) {
      return(NULL)
    }
    
    n<-length(v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)
    if (n!=0) {
      sampleDT <- data.table( 
        fd_method= c(sapply(1:n, function(i) v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods[[i]]@fd_method)),
        fd= c(sapply(1:n, function(i) v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods[[i]]@fd)),
        fd_scale= c(sapply(1:n, function(i) v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods[[i]]@fd_scale)))
      
      print(sampleDT)
      setkey(sampleDT, "fd_method")
      sampleDT <- unique(sampleDT)
      
      name<-NULL
      
      fractalDF<-data.frame(sampleDT[, name:=NULL])
      
      return({xtable(fractalDF)})
    }else{
      return(NULL)
    }
  })
  
  
  
  output$fractalDimensionsFractalUI<-renderUI({
    input$calculateFractalDimensionsButton
    myfractalDimensionsFractalUI()
  })
  myfractalDimensionsFractalUI<-reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@fdAnalysis)) {
      return(list())
    }
    return(h3("Fractal dimensions and scales"))
  })
  
  output$fractalDimensionsFractalPlots_fd2d_isotropic <- renderImage({
    input$calculateFractalDimensionsButton
    fractalDimensionsFractalPlots_fd2d_isotropic()
  }, deleteFile = TRUE)
  output$fractalDimensionsFractalPlots_fd2d_squareincr <- renderImage({
    input$calculateFractalDimensionsButton
    fractalDimensionsFractalPlots_fd2d_squareincr()
  }, deleteFile = TRUE)
  output$fractalDimensionsFractalPlots_fd2d_filter1 <- renderImage({
    input$calculateFractalDimensionsButton
    fractalDimensionsFractalPlots_fd2d_filter1()
  }, deleteFile = TRUE)
  
  
  fractalDimensionsFractalPlots_fd2d_isotropic <- reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@fdAnalysis)) {
      return(list(src = tempfile()))
    }
    if (length(v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)!=0) {
      
      # Create a Progress object
      progressFractal <- shiny::Progress$new()
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progressFractal$close())
      # Create a closure to update progress.
      updateprogressFractal <- function(value = NULL, detail = NULL, message = NULL) {
        if (!is.null(message)) {
          progressFractal$set(message = message, value = 0)
        }else{
          progressFractal$set(value = value, detail = detail)
        }
      }
      updateprogressFractal(message="2/2 - Calculating images", value=0)
      sampleName<-basename(v$AFMImageAnalyser@AFMImage@fullfilename)
      rf2d <- matrix(v$AFMImageAnalyser@AFMImage@data$h, nrow=v$AFMImageAnalyser@AFMImage@samplesperline)
      updateprogressFractal(value= 1/4, detail = "1/4")
      
      outfile1 <- tempfile(fileext='.png')
      print("saving 1")
      png(outfile1, width=400, height=300)
      fd2d_isotropic <- fd.estim.isotropic(rf2d, p.index = 1, direction='hvd+d-', plot.loglog = TRUE, plot.allpoints = TRUE)
      dev.off()
      print("done 1")
      return(list(src = outfile1,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = "fd2d_isotropic")
      )
    }
    return(list(src = tempfile()))
  })
  fractalDimensionsFractalPlots_fd2d_squareincr <- reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@fdAnalysis)) {
      return(list(src = tempfile()))
    }
    if (length(v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)!=0) {
      # Create a Progress object
      progressFractal <- shiny::Progress$new()
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progressFractal$close())
      # Create a closure to update progress.
      updateprogressFractal <- function(value = NULL, detail = NULL, message = NULL) {
        if (!is.null(message)) {
          progressFractal$set(message = message, value = 0)
        }else{
          progressFractal$set(value = value, detail = detail)
        }
      }
      updateprogressFractal(message="2/2 - Calculating images", value=0)
      sampleName<-basename(v$AFMImageAnalyser@AFMImage@fullfilename)
      rf2d <- matrix(v$AFMImageAnalyser@AFMImage@data$h, nrow=v$AFMImageAnalyser@AFMImage@samplesperline)
      updateprogressFractal(value= 2/4, detail = "2/4")
      
      outfile1 <- tempfile(fileext='.png')
      print("saving 2")
      png(outfile1, width=400, height=300)
      fd2d_squareincr <- fd.estim.squareincr(rf2d, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE)
      dev.off()
      print("done 2")
      return(list(src = outfile1,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = "fd2d_squareincr")
      )
    }
    return(list(src = tempfile()))
  })
  fractalDimensionsFractalPlots_fd2d_filter1 <- reactive({
    if (is.null(v$AFMImageAnalyser)||is.null(v$AFMImageAnalyser@fdAnalysis)) {
      return(list(src = tempfile()))
    }
    if (length(v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)!=0) {
      # Create a Progress object
      progressFractal <- shiny::Progress$new()
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progressFractal$close())
      # Create a closure to update progress.
      updateprogressFractal <- function(value = NULL, detail = NULL, message = NULL) {
        if (!is.null(message)) {
          progressFractal$set(message = message, value = 0)
        }else{
          progressFractal$set(value = value, detail = detail)
        }
      }
      updateprogressFractal(message="2/2 - Calculating images", value=0)
      
      sampleName<-basename(v$AFMImageAnalyser@AFMImage@fullfilename)
      rf2d <- matrix(v$AFMImageAnalyser@AFMImage@data$h, nrow=v$AFMImageAnalyser@AFMImage@samplesperline)
      updateprogressFractal(value= 3/4, detail = "3/4")
      
      print("saving 3")
      outfile1 <- tempfile(fileext='.png')
      png(outfile1, width=400, height=300)
      fd2d_filter1 <- fd.estim.filter1(rf2d, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE)
      dev.off()
      print("done 3")
      return(list(src = outfile1,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = "fd2d_filter1")
      )
    }
    return(list(src = tempfile()))
  })
  
  #
  # 
  #
  output$generateCheckReport <- downloadHandler(
    filename = function() { paste(basename(v$AFMImageAnalyser@AFMImage@fullfilename), '-checkReport.pdf', sep='') },
    content = function(file) {
      if(!is.null(v$AFMImageAnalyser)) {
        print("Exporting check report")
        v$AFMImageAnalyser@fullfilename<-file
        generateAFMImageReport(AFMImageAnalyser= v$AFMImageAnalyser, reportFullfilename=file, isCheckReport = TRUE)
        print("done")
      }
    }
  )
  
  output$generateReport <- downloadHandler(
    filename = function() { paste(basename(v$AFMImageAnalyser@AFMImage@fullfilename), '-fullReport.pdf', sep='') },
    content = function(file) {
      if(!is.null(v$AFMImageAnalyser)) {
        print("Exporting full report")
        v$AFMImageAnalyser@fullfilename<-file
        generateAFMImageReport(AFMImageAnalyser= v$AFMImageAnalyser, reportFullfilename=file, isCheckReport = FALSE)
        print("done")
      }
    }
  )
  
  output$alreadyCalculatedPlot <- renderTable({ 
    
    if(is.null(v$AFMImageAnalyser)) {
      shinyjs::disable("generateCheckReport")
      shinyjs::disable("generateReport")
      return()
    }
    
    dataAvailable<-FALSE
    if(!is.null(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms)&&
       length(v$AFMImageAnalyser@variogramAnalysis@directionalVariograms)!=0) {
      checkAvailable<-"Yes"
      dataAvailable<-TRUE
    }
    else checkAvailable<-"No"
    
    if(!is.null(v$AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram)&&
       length(v$AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram)!=0) {
      modelsAvalable<-"Yes"
      dataAvailable<-TRUE
    } 
    else modelsAvalable<-"No"
    
    if(!is.null(v$AFMImageAnalyser@psdAnalysis)&&
       !is.null(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)&&
       length(v$AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)!=0) {
      print(v$AFMImageAnalyser@psdAnalysis)
      psdAvailable<-"Yes"
      dataAvailable<-TRUE
    }
    else psdAvailable<-"No"
    
    if(!is.null(v$AFMImageAnalyser@threeDimensionAnalysis)&&
       !is.null(v$AFMImageAnalyser@threeDimensionAnalysis@f1)&&
       length(v$AFMImageAnalyser@threeDimensionAnalysis@f1)!=0) {
      print(v$AFMImageAnalyser@threeDimensionAnalysis)
      threeDimensionAnalysisAvailable<-"Yes"
      dataAvailable<-TRUE
    }
    else threeDimensionAnalysisAvailable<-"No"
    
    
    if(!is.null(v$AFMImageAnalyser@fdAnalysis)&&
       length(v$AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)!=0) {
      print(v$AFMImageAnalyser@fdAnalysis)
      fdAvailable<-"Yes"
      dataAvailable<-TRUE
    }
    else fdAvailable<-"No"
    
    if (dataAvailable) {
      shinyjs::enable("generateCheckReport")
      shinyjs::enable("generateReport")
    }
    
    availableCalculationsDF<-data.frame(Calculations=c("PSD analysis",
                                                       "Normality and isotropy",
                                                       "Variogram models",
                                                       "Fractal analysis",
                                                       "Networks analysis"), 
                                        available=c(psdAvailable,
                                                    checkAvailable,
                                                    modelsAvalable,
                                                    fdAvailable,
                                                    threeDimensionAnalysisAvailable))
    #print(availableCalculationsDF)
    return({xtable(availableCalculationsDF)})
    
  }, include.rownames=FALSE, include.colnames=TRUE)
  
  #
  # Networks
  #
  
  output$imageNameNetworks<-renderUI({
    imageName<-displayImageName()
    
    if (is.null(imageName)) {
      output$imageNameNetworks<-renderUI(HTML(c("<h4>please select image first</h4>")))
      return(NULL)
    }
    
    heights<-v$AFMImageAnalyser@AFMImage@data$h*input$heightNetworksslider
    heights<-heights+abs(min(heights))
    
    
    print("imageNameNetworks - not null imaheName")
    output$imageNameNetworks<-imageName
    #print(imageName)
  })
  
  
  output$distNetworksPlot <- renderPlot({
    
    if (is.null(v$AFMImageAnalyser)) {
      return()
    }
    
    heights<-v$AFMImageAnalyser@AFMImage@data$h*input$heightNetworksslider
    heights<-heights+abs(min(heights))
    
    bins <- seq(min(heights), max(heights), length.out = 50)
    
    # draw the histogram with the specified number of bins
    hist(heights, breaks = bins, col = 'darkgray', border = 'white')
    updateSliderInput(session, 'filterNetworksslider', value = c(1,10),
                      min = 0.1, max = ceiling(max(heights*110/100)), step = 0.1)
  })
  
  
  output$newImageNetworksPlot <- renderPlot({
    input$checkFilterNetworksButton
    print("checkFilterNetworksButton button pushed")
    
    
    if (is.null(input$checkFilterNetworksButton)) {
      print("input$checkFilterNetworksButton==NULL")
      return(NULL)
    }
    print("input$checkFilterNetworksButton!=NULL")
    if(input$checkFilterNetworksButton == c(0))
    {
      print("input$checkFilterNetworksButton==0")
      return()
    }else{
      
      if (is.null(v$AFMImageAnalyser)) {
        return()
      }
      
      
      
      newAFMImage<-v$AFMImageAnalyser@AFMImage
      
      heights<-newAFMImage@data$h*input$heightNetworksslider
      heights<-heights+abs(min(heights))
      
      heights[heights<input$filterNetworksslider[1]]<-0
      heights[heights>input$filterNetworksslider[2]]<-0
      
      newAFMImage@data$h<-heights
      
      
      
      
      getSpplotFromAFMImage(newAFMImage, expectedWidth=512, expectHeight= 512, withoutLegend=TRUE)
    }
  })
  
  output$skeletonImageNetworksPlot <- renderPlot({
    if (is.null(v$AFMImageAnalyser)||
        is.null(v$AFMImageAnalyser@networksAnalysis)||
        is.null(v$AFMImageAnalyser@networksAnalysis@skeletonGraph)
    ) {
      return()
    }
    library(igraph)
    nbVertices=length(V(v$AFMImageAnalyser@networksAnalysis@skeletonGraph))
    print(nbVertices)
    gridIgraphPlot(AFMImage, v$AFMImageAnalyser@networksAnalysis@skeletonGraph)
    
  })
  
  
  observeEvent(input$calculateNetworksNetworksButton, {
    input$calculateNetworksNetworksButton
    print("calculateNetworksNetworksButton button pushed")
    if (is.null(input$calculateNetworksNetworksButton)) {
      print("input$calculateNetworksNetworksButton==NULL")
      return(NULL)
    }
    print("input$calculateNetworksNetworksButton!=NULL")
    if(input$calculateNetworksNetworksButton == c(0))
    {
      print("input$calculateNetworksNetworksButton==0")
      return()
    }else{
      isolate({
        input$calculateNetworksNetworksButton
        
        # Create a Progress object
        progressCalculateNetworks <- shiny::Progress$new()
        #progressPSD$set(message = "Calculting", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progressCalculateNetworks$close())
        
        print("calculation of Networks")
        
        #createAFMImageAnalyser()
        networkAnalysis<-new ("AFMImageNetworksAnalysis")
        # Create a closure to update progress
        networkAnalysis@updateProgress<- function(value = NULL, detail = NULL, message = NULL) {
          if (!is.null(message)) {
            progressCalculateNetworks$set(message = message, value = 0)
          }else{
            progressCalculateNetworks$set(value = value, detail = detail)
          }
          return(TRUE)
        }
        networkAnalysis@updateProgress(message="Calculating networks", value=0)
        
        
        
        newAFMImage<-copy(v$AFMImageAnalyser@AFMImage)
        
        newAFMImage<-extractAFMImage(newAFMImage,0,0,184)
        
        heights<-newAFMImage@data$h*input$heightNetworksslider
        heights<-heights+abs(min(heights))
        
        heights[heights<input$filterNetworksslider[1]]<-0
        heights[heights>input$filterNetworksslider[2]]<-0
        
        newAFMImage@data$h<-heights
        
        
        networkAnalysis@heightNetworksslider=input$heightNetworksslider
        networkAnalysis@filterNetworkssliderMin=input$filterNetworksslider[1]
        networkAnalysis@filterNetworkssliderMax=input$filterNetworksslider[2]
        
        if(!is.null(v$AFMImageAnalyser@AFMImage)) {
          print("Calculating networks")
          networksAna<-calculateNetworks(AFMImageNetworksAnalysis= networkAnalysis, AFMImage= newAFMImage)
        }
        
        print("calculation of networks done")
        v$AFMImageAnalyser@networksAnalysis<-networksAna
        print("done v$AFMImageAnalyser@networksAnalysis<-networksAna")
      })
    }
    
    
  })
  
  
  
})

