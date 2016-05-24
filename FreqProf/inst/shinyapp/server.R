library(shiny)

shinyServer(function(input, output, session) {

  getDataFromShiny = function(inFile){

    if (is.null(inFile))
      return(NULL)

    # reading a file, whose extension is either csv, bin or fpw,
    # and importing it as a data.frame
    filename = inFile$name

    file.extension = tolower(substr(filename,nchar(filename)-2,nchar(filename)))

    data.behavior = switch(file.extension,
                           csv = read.csv(inFile$datapath),
                           bin = read.bin(inFile$datapath),
                           fpw = read.fpw(inFile$datapath))

    if(is.null(data.behavior)) stop("file extension must be either csv, fpw, or bin")

    # update selected behaviors
    updateCheckboxGroupInput(session, "selected.behaviors",
                             choices = names(data.behavior),
                             selected = input$selected.behaviors)

    data.behavior = data.behavior[,names(data.behavior) %in% input$selected.behaviors]

    if(is.null(ncol(data.behavior))){
      # this means that only one behavior is selected
      dat = as.data.frame(data.behavior)
      names(dat) = input$selected.behaviors
      return(dat)
    }

    if(ncol(data.behavior)>1) return(data.behavior)

    return(NULL)
  }

  getWindowLength = function(unit,window,data){
    if(unit == "bins")
      return(window)
    else
      return(round(window/100*nrow(data)))
  }

  output$distPlot <- renderPlot({
    data.behavior = getDataFromShiny(input$file)
    if(is.null(data.behavior)) return(NULL)

    data.freqprof = freqprof(data.behavior,
                             window = getWindowLength(input$unit_length,input$window,data.behavior),
                             step = input$step,
                             resolution = input$resolution,
                             which = input$which)

    # plotting
    plot_freqprof(data.freqprof,
         gg=input$ggplot,
         panel.in = input$panel.in,
         panel.out = input$panel.out,
         multiPlot = input$multiplot,
         xAxisUnits = input$units,
         tick.every = input$tick.every,
         label.every = input$label.every)
  })

  observe({
    data.behavior = getDataFromShiny(input$file)
    if(is.null(data.behavior)) return(NULL)

    # update range for window length
    if(input$unit_length == "bins"){
      win = round(.25*nrow(data.behavior))
      updateSliderInput(session, "window", value = win,
                        min = 1, max = 4*win, step = 1)
    }
    if(input$unit_length == "percent"){
      updateSliderInput(session, "window", value = 25,
                        min = 1, max = 100, step = 1)
    }

    # update tick.every and label.every
    t.every = round(nrow(data.behavior)/31)
    updateSliderInput(session, "tick.every", value = t.every,
                      min = 1, max = nrow(data.behavior), step = 1)
    updateSliderInput(session, "label.every", value = 3,
                      min = 1, max = 100, step = 1)
  })

  output$downloadData <- downloadHandler(
    filename = "freqprof.csv",
    content = function(file) {
      data.behavior = getDataFromShiny(input$file)
      if(is.null(data.behavior)) return(NULL)

      data.freqprof = freqprof(data.behavior,
                               window = input$window,
                               step = input$step,
                               resolution = input$resolution,
                               which = input$which)

      # which panels will be downloaded?
      panels = c(2)
      if(input$panel.in) panels = c(1,panels)
      if(input$panel.out) panels = c(panels,3)

      write.csv(data.freqprof$data[ data.freqprof$data$panels %in% panels, ], file,row.names=F)
    }
  )

  output$downloadPlotPDF <- downloadHandler(
    filename = function() { paste0("ShinyPlot.pdf") },
    content = function(file) {
      pdf(file,width = input$graphWidth, height = input$graphHeight)
      data.behavior = getDataFromShiny(input$file)
      data.freqprof = freqprof(data.behavior,
                               getWindowLength(input$unit_length,input$window,data.behavior),
                               step = input$step,
                               resolution = input$resolution,
                               which = input$which)
      plot.freqprof(data.freqprof,
                    gg=input$ggplot,
                    panel.in = input$panel.in,
                    panel.out = input$panel.out,
                    multiPlot = input$multiplot,
                    xAxisUnits = input$units)
      dev.off()

      if (file.exists(paste0(file, ".pdf")))
        file.rename(paste0(file, ".pdf"), file)
    })

  output$downloadPlotPNG <- downloadHandler(
    filename = function() { paste0("ShinyPlot.png") },
    content = function(file) {
      png(file,width = input$graphWidth, height = input$graphHeight, units = 'in', res = 100)
      data.behavior = getDataFromShiny(input$file)
      data.freqprof = freqprof(data.behavior,
                               window = getWindowLength(input$unit_length,input$window,data.behavior),
                               step = input$step,
                               resolution = input$resolution,
                               which = input$which)
      plot.freqprof(data.freqprof,
                    gg=input$ggplot,
                    panel.in = input$panel.in,
                    panel.out = input$panel.out,
                    multiPlot = input$multiplot,
                    xAxisUnits = input$units)
      dev.off()

      if (file.exists(paste0(file, ".png")))
        file.rename(paste0(file, ".png"), file)
    })

})


