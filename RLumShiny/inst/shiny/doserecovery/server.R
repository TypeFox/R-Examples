## Server.R
library(Luminescence)
library(shiny)

## read example data set and misapply them for this plot type
data(ExampleData.DeValues, envir = environment())
ExampleData.DeValues <- ExampleData.DeValues$BT998

##############################################################################
###                        MAIN PROGRAM                                    ###
##############################################################################

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
  
  ### GET DATA SETS
  Data<- reactive({
    if(!is.null(datGet())) {
      if(!is.null(datGet2())) {
        data<- list(datGet(), datGet2())
      } else {
        data<- list(datGet())
      }
    } else {
      x.1 <- ExampleData.DeValues[7:11,]
      x.2 <- ExampleData.DeValues[7:11,] * c(runif(5, 0.9, 1.1), 1)
      data<- list(x.1, x.2)
    }
  })
  
  
  output$xlim<- renderUI({
    data<- Data()
    n<- nrow(data[[1]])
    
    sliderInput(inputId = "xlim", label = "Range x-axis", 
                min = 0, max = n*2, 
                value = c(1, n+1))
  })
  
  observe({
    updateTextInput(session, inputId = "xlab", 
                    value = if(input$preheat==TRUE){"Preheat Temperature [\u00B0C]"}else{"# Aliquot"})
  })
  
  #### PLOT ####
  output$main_plot <- renderPlot({
    input$refresh
    
    data<- Data()
    
    
    outputOptions(x = output, name = "xlim", suspendWhenHidden = FALSE)
    validate(
      need(expr = input$xlim, message = 'Waiting for data... Please wait!')
    )
    
    
    # if custom datapoint style get char from separate input panel
    if(input$pch == "custom") {
      pch<- input$custompch
    } else {
      pch<- as.integer(input$pch)-1 #-1 offset in pch values
    }
    # if custom datapoint style get char from separate input panel
    if(input$pch2 == "custom") {
      pch2<- input$custompch2
    } else {
      pch2<- as.integer(input$pch2)-1 #-1 offset in pch values
    }
    
    
    # if custom datapoint color get RGB code from separate input panel
    if(input$color == "custom") {
      color<- input$rgb
    } else {
      color<- input$color
    }
    
    if(length(data) > 1) {
      # if custom datapoint color get RGB code from separate input panel
      if(input$color2 == "custom") {
        color2<- input$rgb2
      } else {
        color2<- input$color2
      }
    } else {
      if(input$preheat == TRUE) {
        color2<- color
      } else {
      color2<- "white"
      }
    }
    
    if(length(data)==1){
      given.dose<- input$dose
      legend<- input$legendname
    } else {
      given.dose<- c(input$dose, input$dose2)
      legend<- c(input$legendname, input$legendname2)
    }
    
    # save all arguments in a list
    args<- list(values = data, 
                error.range = input$error,
                given.dose = given.dose,
                summary = input$stats,
                summary.pos = input$sumpos,
                boxplot = input$boxplot,
                legend = legend,
                legend.pos = input$legend.pos,
                main = input$main,
                mtext = input$mtext,
                col = c(color, color2),
                pch = c(pch, pch2),
                xlab = input$xlab,
                ylab = input$ylab,
                xlim = input$xlim,
                ylim = input$ylim,
                cex = input$cex)
    
    if(input$preheat == TRUE) {
      
      n<- length(data[[1]][,1])
      ph<- c(input$ph1, input$ph2, input$ph3, input$ph4, input$ph5, input$ph6, input$ph7, input$ph8)
      ph<- ph[1:n]
      
      args<- c(args, "preheat" = NA)
      args$preheat<- ph
      
      args$pch<- rep(args$pch, n)
      args$col<- rep(args$col, n)
      
    }
    
    # plot Abanico Plot 
    do.call(what = plot_DRTResults, args = args)
    
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
        
        # plot Abanico Plot 
        do.call(what = plot_DRTResults, args = args)
        
        dev.off()
      },#EO content =,
      contentType = "image"
    )#EndOf::dowmloadHandler()
    
  })
  
  # renderTable() that prints the data to the second tab
  output$dataset<- renderDataTable(
    options = list(pageLength = 10, autoWidth = FALSE),
    callback = "function(table) {
  table.on('click.dt', 'tr', function() {
  $(this).toggleClass('selected');
  Shiny.onInputChange('rows',
  table.rows('.selected').data().toArray());
  });}",
{
  data<- Data()
  colnames(data[[1]])<- c("De", "De error")
  data[[1]]
})##EndOf::renterTable()


# renderTable() that prints the data to the second tab
output$dataset2<- renderDataTable(
  options = list(pageLength = 10, autoWidth = FALSE),
  callback = "function(table) {
  table.on('click.dt', 'tr', function() {
  $(this).toggleClass('selected');
  Shiny.onInputChange('rows',
  table.rows('.selected').data().toArray());
  });}",
{
  data<- Data()
  if(length(data)>1) {
    colnames(data[[2]])<- c("De", "De error")
    data[[2]]
  }
})##EndOf::renterTable()

})
