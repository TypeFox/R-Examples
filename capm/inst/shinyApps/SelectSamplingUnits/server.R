library(rgdal)
shinyServer(function(input, output) {
  
  ## Two-stage cluster sampling.
  
  Universe <- function() {
    if (input$examples) {
      data(psu.ssu)
      return(psu.ssu)
    } else if (is.null(input$psu.ssu)) {
      return()
    } else {
      return(read.csv(input$psu.ssu$datapath,
                      sep = input$sep,
                      quote = input$quote,
                      header = input$header))
    }
  }
  
  ## Title to render above the respective uploaded files.
  FileTitle <- function() {
    if(!is.null(Universe())) {
      return('Uploaded file')
    } else {
      return()
    }
  }
  
  ## Selection of PSUs.
  PSU <- function() {
    if (input$design != 'twostage') {
      return()
    } else {
      if (is.null(Universe()) | !is.numeric(input$psu)) {
        return()
      } else {
        return(SamplePPS(Universe(), input$psu))
      }
    }
  }
  
  ## Selection of SSUs or simple sampling units.
  SSU <- function() {
    if (input$design != 'twostage') {
      return()
    } else {
      if (is.null(PSU()) | !is.numeric(input$ssu)) {
        return()
      } else {
        return(SampleSystematic(psu.ssu = PSU(), su = input$ssu))
      } 
    }
  }
  
  ## Shapefile with the universe of PSUs.
  Maps <- function() {
    if (input$examples) {
      shape.path <- system.file('extdata', package="capm")
      shape.name <- 'santos'
    } else {
      shape.path <- input$shape.path
      shape.name <- input$shape.name
    }
    return(readOGR(shape.path, shape.name))
    
  }
  
  ## Write kml files of the selected PSUs.  
  
  # Redefinition of MapkmlPSU function from capm to write KMLs in a
  # user-defined directory (new argument: write.to.path).
  MapkmlPSU2 <- function (shape = NULL, psu = NULL, id = NULL,
                          path = '.', write.to.path = NULL) 
  {
    if (class(shape) == "SpatialPolygonsDataFrame") {
      tmp <- shape
    }
    else {
      tmp <- readOGR(path, shape)
    }
    tmp <- spTransform(tmp, CRS("+proj=longlat +ellps=WGS84"))
    tmp2 = NULL
    for (i in 1:length(psu)) {
      tmp1 <- tmp[which(as.character(tmp@data[, id]) == psu[i]), 
                  ]
      writeOGR(tmp1, dsn = paste0(write.to.path, '/', eval(psu[i]), ".kml"), 
               layer = 'selected_psu', driver = "KML",
               overwrite_layer = TRUE)
      tmp2[i] <- which(as.character(tmp@data[, id]) == psu[i])
    }
    tmp2 <- tmp[tmp2, ]
    if (file.exists("all_psu.kml")) {
      file.remove("all_psu.kml")
    }
    writeOGR(tmp2, dsn = paste0(write.to.path, '/', "all_psu.kml"),
             layer = 'all_selected_psu', overwrite_layer = TRUE,
             driver = "KML")
  }
  
  # Write the kmls.
  observe({
    if (input$examples) {
      shape.path <- system.file('extdata', package="capm")
      shape.name <- 'santos'
    } else {
      .name <- input$shape.name
    }
    if (input$kml == 0) {
      return()
    }
    isolate({
      if (input$examples) {
        shape.path <- system.file('extdata', package="capm")
      } else {
        shape.path <- input$shape.path
      }
      MapkmlPSU2(shape.name, PSU()[ , 1], 1,
                 path = shape.path, write.to.path = input$write.to.path)
    })
  })
  
  ## Systematic sampling.
  Systematic <- function() {
    if (input$design != 'systematic') {
      return()
    } else {
      if (is.numeric(input$N) & is.numeric(input$su)) {
        return(data.frame('Sampling_units' = SampleSystematic(N = input$N,
                                                              su = input$su)))
      } else {
        return()
      }
    }
  }
  
  ## Stratified sampling.
  Stratified <- function() {
    if (input$design != 'stratified') {
      return()
    } else {
      if (!is.null(input$strata.N) & !is.null(input$strata.su)) {
        N <- as.numeric(strsplit(input$strata.N, ',')[[1]])
        names.N <- strsplit(input$strata.names, ',')[[1]]
        su <- as.numeric(strsplit(input$strata.su, ',')[[1]])
        if (length(N) == length(su) &
              length(N) == length(names.N) & length(N) > 1) {
        names(N) <- names.N
          return(SampleSystematic(N = N, su = su))
        }
      } else {
        return()
      }
    }
  }
    
  output$selected <- renderTable(
    if (!is.null(SSU()) & input$design == 'twostage') {
      SSU()
    } else {
      if (!is.null(Systematic()) & input$design == 'systematic') {
        Systematic()
      } else {
        if (!is.null(Stratified()) & input$design == 'stratified') {
          Stratified()
        }
      }
    },
    digits = 0
  )
  
  output$downloadData <- downloadHandler(
    filename = 'selected_psu.csv',
    content = function(file) {
      write.csv(SSU(), file, row.names = FALSE)
    }
  )
  
  output$file.title <- renderText(FileTitle())
  output$universe <- renderDataTable({
    Universe()
  }, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  
  output$maps <- renderPlot({
    if (input$get.map == 0)
      return()
    
    isolate({
      if (!is.null(input$psu) & !is.null(Universe())) {
        tmp <- NULL
        for (i in 1:length(PSU()[ , 1])) {
          tmp[i] <- which(
            as.character(Maps()@data[, input$col]) == PSU()[ , 1][i])
        }
        plot(Maps(), border = 'grey', axes = T, las = 1,
             xlab = 'Easting', ylab = 'Northing')
        plot(Maps()[tmp, ], col = 'red', border = 'red', add = T)
      }
    })
  })
  
})                    