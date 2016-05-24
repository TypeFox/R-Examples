## Server.R
library(Luminescence)
library(shiny)
library(googleVis)

## MAIN FUNCTION
shinyServer(function(input, output, session) {
  
  observeEvent(input$exit, {
    stopApp(message("Goodbye!"))
  })
  
  # function to convert coordinates to degree decimal format
  coord_conv<- function(x, id) {
    if(id=="degDecMin") {
      x<- paste(as.character(sum(input$degN_1, input$decMinN/60)),":",
                as.character(sum(input$degE_1, input$decMinE/60)), sep = "")
    }
    if(id=="degMinSec") {
      x<- paste(as.character(sum(input$degN_2,input$minN/60,input$secN/3600)),":",
                as.character(sum(input$degE_2,input$minE/60,input$secE/3600)), sep = "")
    }
    return(x)
  }
  
  
  # coordinate conversion
  coords<- reactive({
    if(input$coords != "decDeg") {
      LatLong<- ifelse(input$coords=="degDecMin",
                       coord_conv(, id="degDecMin"),  # YES
                       coord_conv(, id ="degMinSec"))  # NO
    } else {
      LatLong<- paste(input$decDegN,":",input$decDegE,sep="")
    }
    
    # return data frame
    d<- data.frame(LatLong = LatLong, tip = "Site")
    return(d)
  })
  
  # googleVis Map options
  myOptionsMap<- reactive({
    opt<- list(enableScrollWheel = TRUE,
               showTip = TRUE,
               useMapTypeControl = TRUE,
               mapType = "terrain")
    return(opt)
  })
  
  # render googleVis map
  output$map<- renderGvis({
    # refresh plot on button press
    input$refresh
    
    gvisMap(data = coords(), locationvar = "LatLong", tipvar = "tip", options = myOptionsMap())
  })
  
  
  # get results from calc_CosmicDoseRate() function
  get_results<- reactive({
    
    # get coordinates
    coords<- as.vector(coords()$LatLong)
    lat<- as.numeric(unlist(strsplit(x = coords, split = ":"))[1])
    long<- as.numeric(unlist(strsplit(x = coords, split = ":"))[2])
    
    # get absorber properties
    depth<- na.omit(c(input$depth_1, input$depth_2, input$depth_3, input$depth_4, input$depth_5))
    density<- na.omit(c(input$density_1, input$density_2, input$density_3, input$density_4, input$density_5))
    
    t<- get_RLum.Results(calc_CosmicDoseRate(depth = depth, 
                                             density = density, 
                                             latitude = lat, 
                                             longitude = long, 
                                             altitude = input$altitude, 
                                             corr.fieldChanges = input$corr, 
                                             est.age = input$estage, 
                                             half.depth = input$half, 
                                             error = input$error), "summary")
    return(t)
  })
  
  
  # render results for mode 1 and 2
  output$results<- renderUI({
    
    # refresh plot on button press
    input$refresh
    
    if(input$mode == "sAsS" || input$mode == "xAsS") {
      
      t<- get_results()
      HTML(
        if(input$mode=="xAsS") { 
          paste("<font size='2'>","Sample depth: ", "<code>", 
                sum(na.omit(input$depth_1), na.omit(input$depth_2), na.omit(input$depth_3), na.omit(input$depth_4), na.omit(input$depth_5)),
                "m", "</code>", "</font>", "<br>")
        },
        "<font size='2'>","Total absorber: ", "<code>", t$total_absorber.gcm2, "g/cm\u00b2", "</code>", "</font>", "<br>", 
        "<font size='2'>","Cosmic dose rate (uncorrected): ", "<code>", round(t$d0, 3), "Gy/ka", "</code>", "</font>", "<br>", 
        "<font size='2'>","Geomagnetic latitude: ", "<code>", round(t$geom_lat, 2), "\u00b0", "</code>", "</font>", "<br>", 
        "<font size='2'>","Cosmic dose rate (corrected): ", "<code>", round(t$dc, 3),"\u00b1", round(t$dc/100*input$error, 3), "Gy/ka", "</code>", "<br> ", "</font>"
        
      )
    } 
  })
  
  # render results for mode 3
  output$resultsTable<- renderDataTable({
    
    # refresh plot on button press
    input$refresh
    
    if(input$mode == "sAxS") {
      
      t<- get_results()
      
      table<- as.data.frame(cbind(t$depth, t$total_absorber.gcm2, round(t$d0, 3), round(t$dc,3), round(t$dc/100*input$error, 3)))
      colnames(table)<- c("Depth (m)", 
                          "Absorber (g/cm\u00b2)", 
                          "Dc (Gy/ka) [uncorrected]", 
                          "Dc (Gy/ka) [corrected]", 
                          "Dc error (Gy/ka)")
      table
    }
  }, options=list(autoWidth = FALSE, paging = FALSE, processing = TRUE)) # jQuery DataTables options (http://datatables.net)  
})