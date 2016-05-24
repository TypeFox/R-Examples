library(shiny)
library(chipPCR)

# server for the Shiny app
shinyServer(function(input, output) {
  
  processed.data <- reactive({
    dat <- switch(input[["csv.type"]], 
                  csv1 = read.csv(input[["input.file"]][["datapath"]], 
                                  header = input[["header"]]),
                  csv2 = read.csv2(input[["input.file"]][["datapath"]], 
                                   header = input[["header"]]))
    if(!input[["header"]])
      colnames(dat) <- paste0("Column", 1L:ncol(dat))
    dat
  })
  
  
  res.amptest <- reactive({
    dat <- processed.data()
    res <- lapply(1L:ncol(dat), function(i)
      amptester(y = dat[, i], 
                manual = input[["amptester.manual"]]  == "manual",
                noiselevel = input[["amptester.noiselevel"]],
                background = summary(bg.max(1L:length(dat[, i]), 
                                            dat[, i]), print = FALSE)[1:2]))
    res
  })
  
  output[["input.data"]] <- renderTable({
    processed.data()
  })
  
  output[["amptester.summs.plots"]] <- renderUI({
    anal.list <- lapply(1L:length(res.amptest()), function(i) {
      list(plotOutput(paste0("plot", i)), verbatimTextOutput(paste0("summ", i)))
    })
    do.call(tagList, unlist(anal.list, recursive = FALSE))
  })
  
  
  for (i in 1L:300) {
    local({
      my.i <- i
      
      output[[paste0("plot", my.i)]] <- renderPlot(plot(res.amptest()[[my.i]]))
      output[[paste0("summ", my.i)]] <- renderPrint({
        cat(colnames(processed.data())[my.i], "\n")  
        summary(res.amptest()[[my.i]])                                          
      })
    })
  }
  
  
  output[["amptest.summary"]] <- renderUI({
    uiOutput("amptester.summs.plots")
  })
  
  output[["amptest.table"]] <- renderTable({
    dat <- data.frame(t(sapply(res.amptest(), function(i) summary(i, print = FALSE))))
    for (i in c("tht.dec", "slt.dec", "rgt.dec", "shap.noisy", "lrt.test"))
      dat[, i] <- as.logical(dat[, i])
    dat
  })
  
  output[["amptest.data"]] <- renderTable({
    dat <- data.frame(sapply(res.amptest(), function(i) i))
    colnames(dat) <- colnames(processed.data())
    dat
  })
  
  output[["download.table"]] <- downloadHandler(
    filename = "amptester_report.csv",
    content <- function(file) {
      dat <- data.frame(t(sapply(res.amptest(), function(i) summary(i, print = FALSE))))
      for (i in c("tht.dec", "slt.dec", "rgt.dec", "shap.noisy", "lrt.test"))
        dat[, i] <- as.logical(dat[, i])
      #improve it
      write.csv(dat, file)
    }
  )
  
  output[["download.data"]] <- downloadHandler(
    filename = "amptester_data.csv",
    content <- function(file) {
      dat <- data.frame(sapply(res.amptest(), function(i) i))
      colnames(dat) <- colnames(processed.data())
      #improve it
      write.csv(dat, file)
    }
  )
  
  output[["download.result"]] <- downloadHandler(
    filename  = "amptester_report.html",
    content <- function(file) {
      knitr:::knit(input = "amptester_report.Rmd", 
                   output = "amptester_report.md", quiet = TRUE)
      markdown:::markdownToHTML("amptester_report.md", file)
    }
  )
})



