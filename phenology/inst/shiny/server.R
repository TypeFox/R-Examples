library(shiny)
# if (any(installed.packages()[,1]=="shinyIncubator")) library(shinyIncubator)
# runApp(".", launch.browser = TRUE)

library(phenology)

# Define server logic required to plot various variables against mpg

shinyServer(function(input, output, session) {
  

    
    output$resultPlot=renderPlot({
      
      # input$file1 will be NULL initially. After the user selects and uploads a 
      # file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
      # columns. The 'datapath' column will contain the local filenames where the 
      # data can be found.
        
      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      # dest <- "/Users/marc/Documents/Espace_de_travail_R/Phenology/Cayenne/Complete.txt"
      
      dest <- inFile$datapath
      table <- readLines(dest, warn=FALSE)
      
      rp <- .read_phenology(list(table), header=NULL, 
                            reference=NULL, month_ref= NULL, 
                            format=NULL, nm=dest)
      Formated <- add_phenology(previous=NULL, add=rp$DATA, 
                                reference=rp$reference, format=rp$format, silent=TRUE)
      
 #     withProgress(session, min=1, max=5, {setProgress(message = 'Calculation in progress', detail = 'This may take a while...')
        
      pfixed <- c(Min=0, Flat=0)
      Par <- par_init(data=Formated, parametersfixed=pfixed)
 #     if ("package:shinyIncubator" %in% search() ) setProgress(value = 1)
      result1 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
 #     x <- result1$par
 #     pfixed <- c(Flat=0)
 #     Par <- c(x, Min=0.1)
 #     if ("package:shinyIncubator" %in% search() ) setProgress(value = 2)
#      result2 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
#      pfixed <- c(Min=0)
 #     Par <- c(x, Flat=1)
 #     if ("package:shinyIncubator" %in% search() ) setProgress(value = 3)
 #     result3 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
 #     pfixed <- NULL
 #     Par <- c(x, Flat=1, Min=0.1)
#      if ("package:shinyIncubator" %in% search() ) setProgress(value = 4)
#      result4 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
      
#      result <- list(result1, result2, result3, result4)[which.min(c(2*(result1$value+length(result1$par)), 2*(result2$value+length(result2$par)), 2*(result3$value+length(result3$par)), 2*(result4$value+length(result4$par))))][[1]]
#      if ("package:shinyIncubator" %in% search() ) setProgress(value = 5)
        result <- result1
      x <- plot(result, series="all", moon=FALSE, progressbar = FALSE, 
                growlnotify = FALSE)
      
      output$resultsInfo=renderPrint({print(x)})
      
      # End withProgress()
 #     })
     })
    
  
})
