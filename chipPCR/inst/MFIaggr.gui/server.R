library(shiny)
library(chipPCR)


# server for the Shiny app
shinyServer(function(input, output) {
  
  #check if no data is loaded or no example used
  null.input <- reactive({
    is.null(input[["input.file"]]) && input[["run.example"]] == 0
  })
  
  processed.data <- reactive({
    #after loading any file it would be possible to start an example
    if(is.null(input[["input.file"]])) {
      dat <- VIMCFX96_60[, 1L:16]
    } else {
      dat <- switch(input[["csv.type"]], 
                    csv1 = read.csv(input[["input.file"]][["datapath"]], 
                                    header = input[["header"]]),
                    csv2 = read.csv2(input[["input.file"]][["datapath"]], 
                                     header = input[["header"]]))
      if(!input[["header"]])
        colnames(dat) <- paste0("Column", 1L:ncol(dat))
    }
    
    dat
  })
  
  #dabset before and after data input
  output[["dynamic.tabset"]] <- renderUI({
    if(null.input()) {
      tabPanel("No input detected",
               HTML('<p><img src="https://raw.githubusercontent.com/michbur/chipPCR/master/vignettes/logo.png" width="55%" height="55%"/></p>'))
    } else {
      tabsetPanel(
        tabPanel("Input data", tableOutput("input.data")),
        tabPanel("Results with graphics", plotOutput("refMFI.plot"), 
                 verbatimTextOutput("refMFI.summary")),
        tabPanel("Results - table", tableOutput("refMFI.table")),
        tabPanel("All curves plot", plotOutput("allp.plot"))
      )
    }
  })
  
  res.mfi <- reactive({  
    dat <- processed.data()
    
    res <- MFIaggr(x = dat, cyc = input[["cyc.col"]],
                   fluo = (1L:ncol(dat))[-input[["cyc.col"]]], RSD = input[["RSD"]], 
                   rob = input[["rob"]], 
                   llul = c(input[["llul.low"]],input[["llul.up"]]))
    res
  })
  
  output[["input.data"]] <- renderTable({
    processed.data()
  })
  
  output[["refMFI.plot"]] <- renderPlot({
    plot(res.mfi())
  })
  
  output[["allp.plot"]] <- renderPlot({
    dat <- processed.data()
    plotCurves(dat[[input[["cyc.col"]]]], dat[, -input[["cyc.col"]]], CPP = TRUE, type = "l")
  })
  
  
  output[["refMFI.summary"]] <- renderPrint({
    summary(res.mfi())
    
  })
  
  output[["refMFI.table"]] <- renderTable({
    slot(res.mfi(), ".Data")
    
  })
  
  
  output[["download.table"]] <- downloadHandler(
    filename = "refMFI_report.csv",
    content <- function(file) {
      write.csv(slot(res.mfi(), ".Data"), file)
    }
  )
  
  output[["download.result"]] <- downloadHandler(
    filename  = "refMFI_report.html",
    content <- function(file) {
      knitr:::knit(input = "refMFI_report.Rmd", 
                   output = "refMFI_report.md", quiet = TRUE)
      markdown:::markdownToHTML("refMFI_report.md", file)
    }
  )
  
})