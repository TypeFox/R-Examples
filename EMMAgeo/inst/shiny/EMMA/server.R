#' Function to generate GUI server for EMMA.
#' 
#' A graphical user interface (GUI) server backend is started.
#' 
#' @param input A collection of input objects.
#' @param output A collection of input objects.
#' @return A GUI server.
#' @author Michael Dietze, Elisabeth Dietze
#' @keywords EMMA
#' @examples
#' 
#' ## Not run
#' EMMA.GUI()
#' 
#' @export shinyServer

data(X.artificial, envir = environment())
X <- X.artificial

shinyServer(function(
  input, 
  output
  ) {
  
  ## define reactive file loading
  get.data <- reactive({
    
    ## read file
    file.user <- input$file.user
    
    if(is.null(file.user)) {

      return(NULL)
    } else {
     
      return(read.table(file = file.user$datapath, 
                        sep = input$sep)) 
    }
  })
  
  ## render tabular output
  output$tabular <- renderTable({
    
    ## check if file is loaded and overwrite example data
    if(!is.null(get.data())) {
      X <- get.data()
    }
    
    X <- data.frame(x = X)
    
    if(input$Classes == TRUE) {
      classes <- X[1,]
      X <- X[-1,]
    } else {
      classes <- seq(from = 1, to = ncol(X))
    }
    
    if(input$ID == TRUE) {
      ID <- X[,1]
      classes <- classes[-1]
      X <- X[,-1]
    } else {
      ID <- seq(from = 1, to = nrow(X))
    }
    
    if(input$Classes == TRUE) {
      ID <- ID[-1]
    }
    
    X
  })
  
  ## render GSD plot
  output$graphical <- renderPlot({
    
    ## check if file is loaded and overwrite example data
    if(!is.null(get.data())) {
      X <- get.data()
    }
    
    X <- data.frame(x = X)
    
    if(input$Classes == TRUE) {
      classes <- X[1,]
      X <- X[-1,]
    } else {
      classes <- seq(from = 1, to = ncol(X))
    }
    
    if(input$ID == TRUE) {
      ID <- X[,1]
      classes <- classes[-1]
      X <- X[,-1]
    } else {
      ID <- seq(from = 1, to = nrow(X))
    }
    
    if(input$Classes == TRUE) {
      ID <- ID[-1]
    }
    
    if(input$log_x == TRUE) {
      log.data <- "x"
    } else {
      log.data <- ""
    }
    
    ## optionally plot curves
    if(input$show == "curves") {
      
      ## create empty plot
      plot(NA,
           xlim = range(classes), 
           ylim = range(X, na.rm = FALSE),
           log = log.data)
      
      for(i in 1:nrow(X)) {
        lines(x = classes,
              y = X[i,])
      }
      
      ## optionally plot map
    } else if (input$show == "map") {
      fields::image.plot(x = seq(from = 1, to = length(classes)),
                         z = t(as.matrix(X)), 
                         log = log.data, ann = FALSE)
      axis(side = 1)
      axis(side = 2)
      
    } else if (input$show == "EMMA") {
      EMMAgeo::EMMA(X = X,
                    q = input$q,
                    l = input$l,
                    plot = TRUE)
    }
  })

})
