plot_expansion <-
function(exp)
  {
  
     if (class(exp) != "expansion")
    {
        stop(paste(exp, " should be an object of class class 'expansion'.", sep = ""), 
            call. = FALSE)
    }
  
  
  outputN <- exp$NORTH
  outputS <- exp$SOUTH
  outputE <- exp$EAST
  outputW <- exp$WEST
  
    p1 <- gvisLineChart(outputN, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the North", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#006400'}]", 
        backgroundColor = "#D1EEEE"))
    p2 <- gvisLineChart(outputS, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the South", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#0000FF'}]", 
        backgroundColor = "#D1EEEE"))
    p3 <- gvisLineChart(outputE, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the East", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#8B0000'}]", 
        backgroundColor = "#D1EEEE"))
    p4 <- gvisLineChart(outputW, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the West", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#FF4500'}]", 
        backgroundColor = "#D1EEEE"))
    ln1 <- gvisMerge(p1, p2, horizontal = TRUE)
    ln2 <- gvisMerge(p3, p4, horizontal = TRUE)
    ln.final <- gvisMerge(ln1, ln2, horizontal = FALSE)
    ln.final$html$footer <- "\n<!-- htmlFooter -->\n<span> \n  R version 3.0.1 (2013-05-16) &#8226; <a href=\"http://code.google.com/p/google-motion-charts-with-r/\">googleVis-0.4.7</a>\n  &#8226; LandSim-0.1.0\n  &#8226; <a href=\"https://developers.google.com/terms/\">Google Terms of Use</a> &#8226; <a href=\"https://google-developers.appspot.com/chart/interactive/docs/gallery/linechart.html#Data_Policy\">Data Policy</a>\n</span></div>\n</body>\n</html>\n"
    output <- list(NORTH = outputN, SOUTH = outputS, EAST = outputE, WEST = outputW)
    plot(ln.final)


  }
