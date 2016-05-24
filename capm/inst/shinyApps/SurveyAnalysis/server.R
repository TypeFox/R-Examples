library(survey)
shinyServer(function(input, output) {
  
  Universe <- function() {
    if (input$design != 'twostage') {
      return()
    } else if (input$examples.twostage) {
      data(psu.ssu)
      return(psu.ssu)
    } else if (is.null(input$psu.ssu)) {
      return()
    } else {
      return(read.csv(input$universe$datapath,
                      sep=input$sep.two1,
                      quote=input$quote.two1,
                      header = input$header))
    }
  }
  
  SampleData <- function() {
    if (input$design == 'twostage') {
      if (input$examples.twostage) {
        data(survey.data)
        return(sapply(survey.data, as.character))
      } else if (is.null(input$sample)) {
        return()
      } else {
        return(sapply(read.csv(input$sample$datapath,
                        sep=input$sep.two2,
                        quote=input$quote.two2,
                        header = input$header), as.character)) 
      }
    } else if (input$design == 'systematic') {
      if (input$example.systematic) {
        data(survey.data)
        return(survey.data[ , -c(1:2)])
      } else if (is.null(input$sample)) {
        return()
      } else {
        return(sapply(read.csv(input$sample$datapath,
                        sep=input$sep.syst,
                        quote=input$quote.syst,
                        header = input$header), as.character))
      }
    } else if (input$design == 'stratified') {
      if (input$example.stratified) {
        data(survey.data)
        strat <- survey.data[ , -c(1:2)]
        strat$strat <- 'Urban'
        strat$strat[round(runif(5, 1, nrow(strat)))] <- 'Rural'
        strat$strat.size <- 144000
        strat$strat.size[strat$strat == 'Rural'] <- 600
        return(strat)
      } else if (is.null(input$sample)) {
        return()
      } else {
        return(sapply(read.csv(input$sample$datapath,
                        sep=input$sep.strat,
                        quote=input$quote.strat,
                        header = input$header), as.character))
      }
    }
  }
  
  FileTitle <- function() {
    if(is.null(Universe()) & is.null(SampleData())) {
      return() 
    } else if (!is.null(Universe()) & !is.null(SampleData())) {
      return(list('Universe:', 'Survey data:'))
    } else if (is.null(Universe()) & !is.null(SampleData())) {
      return(list(NULL, 'Survey data'))
    }
  }
  
  Pyramid <- function() {
    if (is.null(SampleData()) | is.na(input$age.col) | is.na(input$sex.col)) {
      return()
    } else if (is.na(input$cas.col)) {
      return(PlotPopPyramid(SampleData(), input$age.col,
                            input$sex.col))
    } else {
      return(PlotPopPyramid(SampleData(), input$age.col,
                            input$sex.col, input$cas.col))
    }
  }
  
  Design <- function() {
    if (input$design == 'twostage') {
      if (is.null(Universe()) | is.null(SampleData()) |
            is.null(input$psu.col) | is.null(input$ssu.col) |
            is.null(input$psu.2cdl)) {
        return()
      } else {
        return(DesignSurvey(sample = SampleData(), psu.ssu = Universe(),
                            psu.col = input$psu.col, ssu.col = input$ssu.col,
                            psu.2cd = input$psu.2cdl))
      }
    } else if (input$design == 'systematic') {
      if (is.null(SampleData()) | is.na(input$N)) {
        return()
      } else {
        return(DesignSurvey(sample = SampleData(), N = input$N)) 
      }
    } else if (input$design == 'stratified') {
      if (is.null(SampleData()) | input$strata.sizes == '' |
            input$strata.membership == '') {
        return()
      } else {
        strata.sizes <- input$strata.sizes
        strata.membership <- input$strata.membership
        if (!is.na(as.numeric(strata.sizes))) {
          strata.sizes <- as.numeric(strata.sizes)
        }
        if (!is.na(as.numeric(strata.membership))) {
          strata.membership <- as.numeric(strata.membership)
        }
        return(DesignSurvey(sample = SampleData(), N = strata.sizes,
                              strata = strata.membership))
        }
      }
  }
  
  Variables <- function() {
    if (is.null(input$variables)) {
      return()
    } else {
      no.spaces <- gsub(' +', '', input$variables)
      return(unlist(strsplit(no.spaces, ',')))
    }
  }
  
  output$universe <- renderDataTable({
    Universe()
  }, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  
  output$sample <- renderDataTable({
    SampleData()
  }, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  
  output$file.title1 <- renderText(FileTitle()[[1]])
  output$file.title2 <- renderText(FileTitle()[[2]])
  
  output$summ.sample <- renderPrint({
    if (!is.null(SampleData()) & !is.na(input$psu.col) &
          !is.na(input$ssu.col)) {
      summary(SampleData()[ , -c(1:2)])
    }
  })
  
  output$pyramid <- renderPlot({
    Pyramid()
  })
  
  output$variables <- renderTable({
    if (is.null(Design()) | is.null(Variables())) {
      return()
    } else {
      des.var <- names(Design()$variables)
      cbind(Variables = des.var, 'Tipe of estimate' = Variables())
    }
  })
  
  output$estimates <- renderTable({
    if (input$calc.estimates == 0) {
      return()
    }
    isolate({
      if (is.null(Design()) | is.null(Variables())) {
        return()
      } else {
        SummarySurvey(Design(), Variables(), conf.level = input$conf.level)
      }
    })
  })
  
})