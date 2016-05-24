## Server.R
library(Luminescence)
library(shiny)

# load example data
data(ExampleData.DeValues)
data <- ExampleData.DeValues$CA1

## MAIN FUNCTION
shinyServer(function(input, output, session) {
  
  observeEvent(input$exit, {
    stopApp(message("Goodbye!"))
  })
  
  # check and read in file (DATA SET 1)
  datGet<- reactive({
    inFile<- input$file1
    if(is.null(inFile)) return(NULL) # if no file was uploaded return NULL
    return(read.table(file = inFile$datapath, # inFile[1] contains filepath 
                      sep = input$sep, 
                      quote = "", 
                      header = input$headers)) # else return file
  })
  
  
  # dynamically inject sliderInput for x-axis range
  output$xlim<- renderUI({
    
    # check if file is loaded
    # # case 1: yes -> slinderInput with custom values
    if(!is.null(datGet())) {
      
      data<- datGet()
      xlim.plot<- range(hist(data[,1], plot = FALSE)$breaks)
      
      sliderInput(inputId = "xlim", 
                  label = "Range x-axis",
                  min = xlim.plot[1]*0.5, 
                  max = xlim.plot[2]*1.5,
                  value = c(xlim.plot[1], xlim.plot[2]), round=FALSE, step=0.0001)
    }
    
    else { #case 2: no -> sliderInput for example data
      
      xlim.plot<- range(hist(data[,1], plot = FALSE)$breaks)
      
      sliderInput(inputId = "xlim", 
                  label = "Range x-axis",
                  min = xlim.plot[1]*0.5, 
                  max = xlim.plot[2]*1.5,
                  value = c(xlim.plot[1], xlim.plot[2]), round=FALSE, step=0.0001)
    }
  })## EndOf::renderUI()
  
  
  output$main_plot <- renderPlot({
    
    # refresh plot on button press
    input$refresh
    
    # progress bar
    progress<- Progress$new(session, min = 0, max = 3)
    progress$set(message = "Calculation in progress",
                 detail = "Retrieve data")
    on.exit(progress$close())
    
    # make sure that input panels are registered on non-active tabs.
    # by default tabs are suspended and input variables are hence
    # not available
    outputOptions(x = output, name = "xlim", suspendWhenHidden = FALSE)
    
    # check if file is loaded and overwrite example data
    if(!is.null(datGet())) {
      data<- datGet()
    }
    
    progress$set(value = 1)
    progress$set(message = "Calculation in progress",
                 detail = "Get values")
    
    # check if any summary stats are activated, else NA
    if (input$summary) {
      summary<- input$stats
    } else {
      summary<- NA
    }
    
    # if custom datapoint color get RGB code from separate input panel
    if(input$pchColor == "custom") {
      pch.color<- input$pchRgb
    } else {
      pch.color<- input$pchColor
    }
    
    # if custom datapoint color get RGB code from separate input panel
    if(input$barsColor == "custom") {
      bars.color<-  adjustcolor(col = input$barsRgb,
                                alpha.f = input$alpha.bars/100)
    } else {
      bars.color<-  adjustcolor(col = input$barsColor,
                                alpha.f = input$alpha.bars/100)
    }
    
    # if custom datapoint color get RGB code from separate input panel
    if(input$rugsColor == "custom") {
      rugs.color<- input$rugsRgb
    } else {
      rugs.color<- input$rugsColor
    }
    
    # if custom datapoint color get RGB code from separate input panel
    if(input$normalColor == "custom") {
      normal.color<- input$normalRgb
    } else {
      normal.color<- input$normalColor
    }
    
    # update progress bar
    progress$set(value = 2)
    progress$set(message = "Calculation in progress",
                 detail = "Combine values")
    
    colors<- c(bars.color, rugs.color, normal.color, pch.color)
    
    # if custom datapoint style get char from separate input panel
    if(input$pch == "custom") {
      pch<- input$custompch
    } else {
      pch<- as.integer(input$pch)-1 #-1 offset in pch values
    }
    
    # validate(need()) makes sure that all data are available to
    # renderUI({}) before plotting and will wait until there
    validate(
      need(expr = input$xlim, message = 'Waiting for data... Please wait!')
    )
    
    progress$set(value = 3)
    progress$set(message = "Calculation in progress",
                 detail = "Ready to plot")
    
    plot_Histogram(data = data,
                   na.exclude = input$naExclude, 
                   cex.global = input$cex, 
                   pch = pch,
                   xlim = input$xlim,
                   summary.pos = input$sumpos, 
                   mtext = input$mtext, 
                   main = input$main,
                   rug = input$rugs, 
                   se = input$errorBars, 
                   normal_curve = input$norm, 
                   summary = summary,
                   xlab = input$xlab,
                   ylab = c(input$ylab1, input$ylab2),
                   colour = colors)
    
    
    # char vector for code output
    verb.summary<- "c('"
    for(i in 1:length(summary)){
      verb.summary<- paste(verb.summary, summary[i], "','", sep="")
      if(i == length(summary)) {
        verb.summary<- substr(verb.summary, 1, nchar(verb.summary)-2)
        verb.summary<- paste(verb.summary, ")", sep="")
      }
    }
    
    # char vectors for code output
    str1 <- paste("plot_Histogram(data = data, ", sep = "")          
    str2 <- paste("summary.pos = '", input$sumpos,"',", sep = "")      
    str3 <- paste("summary = ",verb.summary,",", sep = "")            
    str4 <- paste("colour = c('",colors[1],"','",colors[2],"','",colors[3],"','",colors[4],"'),", sep = "") 
    str5 <- paste("pch = ",pch,",", sep = "")                         
    str6<- paste("normal_curve = ",input$norm,",", sep = "")            
    str7<- paste("se = ", input$errorBars, ",", sep = "")                      
    str8<- paste("rug = ", input$rugs, ",", sep = "")                   
    str9 <- paste("main = '",input$main,"',", sep = "")               
    str10 <- paste("cex.global = ", input$cex, ",", sep = "")          
    str11 <- paste("mtext = '",input$mtext,"',", sep = "")             
    str12 <- paste("na.exclude = ",input$naExclude,",", sep = "")     
    str13 <- paste("xlab = '", input$xlab,"',", sep="")                
    str14 <- paste("ylab = c('",input$ylab1,"','", input$ylab2,"'),",sep ="") 
    str15 <- paste("xlim = c(", input$xlim[1],",",input$xlim[2],"))", sep="") 
    
    if(input$sep == "\t") { verb.sep<-  "\\t"}
    else {
      verb.sep<- input$sep
    }
    
    str0.1 <- paste("data <- read.delim(file, header = ",input$headers, ", sep= '", verb.sep,"')",
                    sep = "")
    
    str0 <- paste("# To reproduce the plot in your local R environment",
                  "# copy and run the following code to your R console.",
                  "library(Luminescence)",
                  "file<- file.choose()",
                  str0.1,
                  "\n",
                  str1,
                  sep = "\n")
    
    code.output<- paste(str0,
                        str2, str3, str4, str5, str6, str7, str8, str9, str10, 
                        str11, str12, str13, str14, str15,
                        sep="\n   ")
    
    # nested renderText({}) for code output on "R plot code" tab
    output$plotCode<- renderText({
      
      code.output
      
    })##EndOf::renderText({})
    
    output$exportScript <- downloadHandler(
      filename = function() { paste(input$filename, ".", "R", sep="") },
      content = function(file) {
        write(code.output, file)
      },#EO content =,
      contentType = "text"
    )#EndOf::dowmloadHandler()
    
    
    # nested downloadHandler() to print plot to file
    output$exportFile <- downloadHandler(
      filename = function() { paste(input$filename, ".", input$fileformat, sep="") },
      content = function(file) {
        
        # determine desired fileformat and set arguments
        if(input$fileformat == "pdf") {
          pdf(file, 
              width = input$imgwidth, 
              height = input$imgheight, 
              paper = "special",
              useDingbats = FALSE, 
              family = input$fontfamily)
        }
        if(input$fileformat == "svg") {
          svg(file, 
              width = input$imgwidth, 
              height = input$imgheight, 
              family = input$fontfamily)
        }
        if(input$fileformat == "eps") {
          postscript(file, 
                     width = input$imgwidth, 
                     height = input$imgheight, 
                     paper = "special", 
                     family = input$fontfamily)
        }
        
        plot_Histogram(data = data,
                       na.exclude = input$naExclude, 
                       cex.global = input$cex, 
                       pch = pch,
                       xlim = input$xlim,
                       summary.pos = input$sumpos, 
                       mtext = input$mtext, 
                       main = input$main,
                       rug = input$rugs, 
                       se = input$errorBars, 
                       normal_curve = input$norm, 
                       summary = summary,
                       xlab = input$xlab,
                       ylab = c(input$ylab1, input$ylab2),
                       colour = colors)
        
        dev.off()
        
      },#EO content =,
      contentType = "image"
    )#EndOf::dowmloadHandler()
  })##EndOf::renderPlot({})
  
  
  # renderTable() that prints the data to the second tab
  output$dataset<- renderDataTable(
    options = list(pageLength = 10, autoWidth = FALSE),
    callback = "function(table) {
    table.on('click.dt', 'tr', function() {
    $(this).toggleClass('selected');
    Shiny.onInputChange('rows',
    table.rows('.selected').data().toArray());
    });
}",
{
  if(!is.null(datGet())) {
    data<- datGet()
    colnames(data)<- c("De","De error")
    data
  } else {
    colnames(data)<- c("De","De error")
    data
  }
  })##EndOf::renterTable()
  
  
  # reactive function for gVis plots that allow for dynamic input!
  myOptionsCAM<- reactive({
    options<- list(
      page="enable",
      width="500px",
      sort="disable")
    return(options)
  })
  
  # renderTable() to print the results of the
  # central age model (CAM)
  output$CAM<- renderDataTable(
    options = list(pageLength = 10, autoWidth = FALSE),
{
  if(!is.null(datGet())) {
      data<- list(datGet())
  } else {
    data<- list(data)
  }
  t<- as.data.frame(matrix(nrow = length(data), ncol = 7))
  colnames(t)<- c("Data set","n", "log data", "Central dose", "SE abs.", "OD (%)", "OD error (%)")
  res<- lapply(data, function(x) { calc_CentralDose(x, verbose = FALSE, plot = FALSE) })
  for(i in 1:length(res)) {
    t[i,1]<- ifelse(i==1,"pimary","secondary")
    t[i,2]<- length(res[[i]]@data$data[,1])
    t[i,3]<- res[[i]]@data$args$log
    t[i,4:7]<- round(res[[i]]@data$summary[1:4],2)
  }
  t
})##EndOf::renterTable()

})##EndOf::shinyServer(function(input, output)