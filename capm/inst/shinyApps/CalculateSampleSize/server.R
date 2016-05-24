shinyServer(function(input, output) {
  
  # Two-stage designs
  Universe <- function() {
    if (input$design != 'twostage') {
      return()
    } else {
      if (input$examples) {
        data(psu.ssu)
        return(psu.ssu)
      } else if (is.null(input$psu.ssu)) {
        return()
      } else {
        return(read.csv(input$psu.ssu$datapath,
                        sep=input$sep,
                        quote=input$quote,
                        header = input$header))
      }
    }
  }
  
  SampleData <- function() {
    if (input$examples) {
      data(pilot)
      return(pilot)
    } else if (is.null(input$psu.x)) {
      return()
    } else {
      return(read.csv(input$psu.x$datapath,
                      sep=input$sep,
                      quote=input$quote,
                      header = input$header))
    }
  }
  
  FileTitle <- function() {
    if(!is.null(Universe()) & !is.null(SampleData())) {
      return(c('First uploaded file:', 'Second uploaded file:')) 
    } else {
      return()
    }
  }
  
  # Simple design
  Systematic <- function() {
    if(input$design != 'systematic') {
      return()
    } else{
      if(is.numeric(input$N) & is.numeric(input$expected.mean) &
           is.numeric(input$expected.var)) {
        ((input$N - 1) / input$N * input$expected.var) / input$expected.mean^2
      } else {
        return()
      }
    }
  }
  
  # Stratified design
  Stratified <- function() {
    if (input$design != 'stratified') {
      return()
    } else {
      if (input$strata.names != '' & input$strata.N != '' &
            input$strata.mean  != '' & input$strata.var != '') {
        N <- as.numeric(strsplit(input$strata.N, ',')[[1]])
        names(N) <- strsplit(input$strata.names, ',')[[1]]
        expected.mean <- as.numeric(strsplit(input$strata.mean, ',')[[1]])
        expected.var <- as.numeric(strsplit(input$strata.var, ',')[[1]])
        return(cbind(N, expected.mean, expected.var))
      } else{
        return()
      }
    }
  }
  
  output$size <- renderTable({
    if (!is.null(Universe()) & !is.null(SampleData())) {
      Calculate2StageSampleSize(Universe(), SampleData(), input$level,
                                input$error, input$cost, input$min.ssu)
    } else {
      if (!is.null(Systematic())) {
        matrix(c(CalculateSimpleSampleSize(Systematic(), input$N,
                                           input$level,
                                           input$error)),
               dimnames = list('Sample size:', 'Value'))
      } else {
        if (!is.null(Stratified())) {
          CalculateStratifiedSampleSize(Stratified(),
                                        conf.level = input$level,
                                        error = input$error)
        } else {
          return() 
        }
      }
    }
  }, align = c('l', 'r'))
  
  output$file.title1 <- renderText(FileTitle()[1])
  output$universe <- renderDataTable({
    Universe()
  }, options = list(aLengthMenu = c(5, 30, 50), iDisplayLength = 5))
  
  output$file.title2 <- renderText(FileTitle()[2])
  output$sample.data <- renderDataTable({
    SampleData()
  }, options = list(aLengthMenu = c(5, 30, 50), iDisplayLength = 5))
  
})              