## Server.R
library(Luminescence)
library(shiny)

# load example data
data(ExampleData.DeValues)
data <- list(ExampleData.DeValues$CA1, ExampleData.DeValues$CA1)

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
  
  # check and read in file (DATA SET 2)
  datGet2<- reactive({
    inFile2<- input$file2
    if(is.null(inFile2)) return(NULL) # if no file was uploaded return NULL
    return(read.table(file = inFile2$datapath, # inFile[1] contains filepath 
                      sep = input$sep, 
                      quote = "", 
                      header = input$headers)) # else return file
  })
  
  
  # dynamically inject sliderInput for x-axis range
  output$xlim<- renderUI({
    
    # check if file is loaded
    # # case 1: yes -> slinderInput with custom values
    if(!is.null(datGet())) {
      if(!is.null(datGet2())) {
        
        data<- rbind(datGet(),datGet2())
        
        sliderInput(inputId = "xlim", 
                    label = "Range x-axis",
                    min = min(data[,1])*0.25,
                    max = max(data[,1])*1.75,
                    value = c(min(data[,1])*0.9, max(data[,1])*1.1))
        
      } else {
        data<- datGet()
        
        sliderInput(inputId = "xlim", 
                    label = "Range x-axis",
                    min = min(data[,1])*0.25,
                    max = max(data[,1])*1.75,
                    value = c(min(data[,1])*0.9, max(data[,1])*1.1))
      }
    }
    
    else { #case 2: no -> sliderInput for example data
      
      sliderInput(inputId = "xlim", 
                  label = "Range x-axis",
                  min = min(data[[1]][,1])*0.25, 
                  max = max(data[[1]][,1])*1.75,
                  value = c(min(data[[1]][,1])*0.9, max(data[[1]][,1]))*1.05,
                  step = 1, round = 0)
    }
  })## EndOf::renderUI()
  
  # dynamically inject sliderInput for KDE bandwidth
  output$bw<- renderUI({
    
    # check if file is loaded
    # # case 1: yes -> slinderInput with custom values
    if(!is.null(datGet())) {
      if(!is.null(datGet2())) {
        data<- rbind(datGet(),datGet2())
        
        sliderInput(inputId = "bw", 
                    label = "KDE bandwidth", 
                    min = bw.nrd0(data[,1])/4, 
                    max = bw.nrd0(data[,1])*4,
                    value = bw.nrd0(data[,1]))
        
      } else {
        
        data<- datGet()
        
        sliderInput(inputId = "bw", 
                    label = "KDE bandwidth", 
                    min = bw.nrd0(data[,1])/4, 
                    max = bw.nrd0(data[,1])*4,
                    value = bw.nrd0(data[,1]))
      }
    }
    
    else { #case 2: no -> sliderInput for example data
      # logged data
      
      sliderInput(inputId = "bw", 
                  label = "KDE bandwidth", 
                  min = 1, max = 100,
                  value = 15, step = 1)
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
    outputOptions(x = output, name = "bw", suspendWhenHidden = FALSE)
    
    # check if file is loaded and overwrite example data
    if(!is.null(datGet())) {
      data<- list(datGet(), datGet())
    }
    if(!is.null(datGet2())) {
      data2<- datGet2()
    }
    
    if(is.null(datGet()) == FALSE && is.null(datGet2()) == FALSE) {
      data<- datGet()
      data2<- datGet2()
      data<- list(data, data2)
    }
    
    progress$set(value = 1)
    progress$set(message = "Calculation in progress",
                 detail = "Get values")
    
    
    # check if any summary stats are activated, else NA
    if (input$summary) {
      summary<- input$stats
    } else {
      summary<- ""
    }
    
    if(input$logx == TRUE) {
      logx<- "x"
    } else {
      logx<- ""
    }
    
    # if custom polygon color get RGB from separate input panel or "none"
    if(input$polygon == "custom") {
      polygon.col<- adjustcolor(col = input$rgbPolygon, 
                                alpha.f = input$alpha.polygon/100)
    } else {
      if(input$polygon == "none") {
        polygon.col<- "white"
      } else {
        polygon.col<- adjustcolor(col = input$polygon, 
                                  alpha.f = input$alpha.polygon/100) 
      }
    }
    
    # if custom polygon color get RGB from separate input panel or "none"
    if(input$polygon2 == "custom") {
      polygon.col2<- adjustcolor(col = input$rgbPolygon2, 
                                 alpha.f = input$alpha.polygon/100)
    } else {
      if(input$polygon2 == "none") {
        polygon.col2<- "white"
      } else {
        polygon.col2<- adjustcolor(col = input$polygon2, 
                                   alpha.f = input$alpha.polygon/100) 
      }
    }
    
    # update progress bar
    progress$set(value = 2)
    progress$set(message = "Calculation in progress",
                 detail = "Combine values")
    
    poly.col<- c(polygon.col, polygon.col2)
    
    # if custom datapoint color get RGB code from separate input panel
    if(input$color == "custom") {
      color<- input$rgb
    } else {
      color<- input$color
    }
    
    if(!is.null(datGet2())) {
      # if custom datapoint color get RGB code from separate input panel
      if(input$color2 == "custom") {
        color2<- input$rgb2
      } else {
        color2<- input$color2
      }
    } else {
      color2<- adjustcolor("white", alpha.f = 0)
    }
    
    
    # validate(need()) makes sure that all data are available to
    # renderUI({}) before plotting and will wait until there
    validate(
      need(expr = input$xlim, message = 'Waiting for data... Please wait!'),
      need(expr = input$bw, message = 'Waiting for data... Please wait!')
    )
    
    progress$set(value = 3)
    progress$set(message = "Calculation in progress",
                 detail = "Ready to plot")    
    
    plot_KDE(data = data, 
             centrality = input$centrality, 
             cex = input$cex, 
             log = logx,
             xlab = input$xlab,
             ylab = c(input$ylab1, input$ylab2),
             main = input$main,
             weights = input$weights,
             values.cumulative = input$cumulative,
             na.exclude = input$naExclude, 
             dispersion = input$dispersion, 
             summary = summary,
             summary.pos = input$sumpos,
             bw = input$bw,
             xlim = input$xlim,
             polygon.col = poly.col,
             col = c(color, color2))
    
    
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
    str1 <- paste("plot_KDE(data = data, ", sep = "")
    str2 <- paste("centrality = '",input$centrality,"',", sep = "")
    str3 <- paste("bw = ",input$bw,",", sep = "")
    str4 <- paste("dispersion = '",input$dispersion,"',", sep = "")
    str5 <- paste("log = '",logx, "',", sep = "")
    str6 <- paste("col = c('",color2,"','",color,"'),", sep = "")
    str7 <- paste("main = '",input$main,"',", sep = "")
    str8 <- paste("cex = ", input$cex, ",", sep = "")
    str9 <- paste("summary = ",verb.summary,",", sep = "")
    str10 <- paste("polygon.col = c('",polygon.col,"','",polygon.col2,"'),", sep = "")
    str11 <- paste("na.exclude = ",input$naExclude,",", sep = "")
    str12 <- paste("ylab = c('", input$ylab1,"','",input$ylab2,"'),", sep="")
    str13 <- paste("xlab = '",input$xlab,"',",sep ="")
    str14 <- paste("xlim = c(", input$xlim[1],",",input$xlim[2],"),", sep="")
    str15 <- paste("weights = ", input$weights, ",", sep = "")
    str16 <- paste("values.cumulative = ", input$cumulative, ",", sep = "")
    str17 <- paste("summary.pos = '", input$sumpos, "')", sep = "")
    
    verb.sep<- ifelse(input$sep == "\t", "\\t", ifelse(is.null(input$sep), "\\t", input$sep))
    
    str0.1 <- paste("data <- read.delim(file, header = ",input$headers, ", sep= '", verb.sep,"')",
                    sep = "")
    if(!is.null(datGet2())) {
      str0.2.0 <- "file2<- file.choose()"
      str0.2.1 <- paste("data2 <- read.delim(file2, header = ",input$headers, ", sep= '", verb.sep,"')",
                        sep= "")
      str0.2.2 <- "data<- list(data, data2)"
      str0.1 <- paste(str0.1, str0.2.0, str0.2.1, str0.2.2, sep = "\n")
    }
    
    str0 <- paste("# To reproduce the plot in your local R environment",
                  "# copy and run the following code to your R console.",
                  "library(Luminescence)",
                  "file<- file.choose()",
                  str0.1,
                  "\n",
                  str1,
                  sep = "\n")
    
    code.output<- paste(str0, str2, str3, str4, str5, 
                        str6, str7, str8, str9, str10, str11, str12,
                        str13, str14, str15, str16, str17,
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
        
        if(is.null(input$fileformat))       
          updateRadioButtons(session, inputId = "fileformat", label = "Fileformat", selected = "pdf")
        
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
        
        plot_KDE(data = data, 
                 centrality = input$centrality, 
                 cex = input$cex, 
                 log = logx,
                 xlab = input$xlab,
                 ylab = c(input$ylab1, input$ylab2),
                 main = input$main,
                 weights = input$weights,
                 values.cumulative = input$cumulative,
                 na.exclude = input$naExclude, 
                 dispersion = input$dispersion, 
                 summary = summary,
                 summary.pos = input$sumpos,
                 bw = input$bw,
                 xlim = input$xlim,
                 polygon.col = poly.col,
                 col = c(color2, color))
        
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
    data<- data[[1]]
    colnames(data)<- c("De","De error")
    data
  }
})##EndOf::renterTable()

# renderTable() that prints the secondary data to the second tab
output$dataset2<- renderDataTable(
  options = list(pageLength = 10, autoWidth = FALSE),
  callback = "function(table) {
  table.on('click.dt', 'tr', function() {
  $(this).toggleClass('selected');
  Shiny.onInputChange('rows',
  table.rows('.selected').data().toArray());
  });
  }",
{
  if(!is.null(datGet2())) {
    data<- datGet2()
    colnames(data)<- c("De","De error")
    data
  } else {
  }
})##EndOf::renterTable()


# renderTable() to print the results of the
# central age model (CAM)
output$CAM<- renderDataTable(
  options = list(pageLength = 10, autoWidth = FALSE),
{
  if(!is.null(datGet())) {
    if(!is.null(datGet2())) {
      data<- list(datGet(), datGet2())
    } else {
      data<- list(datGet())
    }
  } else {
    data<- list(data[[1]])
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