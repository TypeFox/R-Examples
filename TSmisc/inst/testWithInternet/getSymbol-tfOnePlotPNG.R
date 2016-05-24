require("TSmisc")  
require("tfplot")  

fred <- TSconnect("getSymbol", dbname="FRED") 

tmp <- tempdir()

files <- paste(tmp,
  c("US-M2small.png", "US-CPIsmall.png","Fordsmall.png"),sep="/")

png(file=files[1], width=480, height=240, pointsize=12, bg = "white")
# mv US-M2small.png US-M2small.png.orig ; pngcrush US-M2small.png.orig US-M2small.png
#png(file="US-M2.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(percentChange(TSget(serIDs="M2", con=fred), lag=52), 
    Title = "Running commentary, blah, blah, blah", 
    subtitle="Broad Money (M2)",
    ylab= "y/y percent change*",
    source="Source: Federal Reserve Bank of St.Louis (M2)",
    footnoteLeft = "seasonally adjusted data",
    footnoteRight = "* approximated by 52 week growth",
    lastObs = TRUE )
dev.off()

png(file=files[2],width=480, height=240, pointsize=12, bg = "white")
# mv US-CPIsmall.png US-CPIsmall.png.orig ; pngcrush US-CPIsmall.png.orig US-CPIsmall.png
#png(file="US-CPI.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(annualizedGrowth(TSget("CPIAUCNS", fred)), start=c(2005,1),
    Title = "U.S. Consumer Price Growth, m/m at Annual Rates", 
    subtitle="Consumer Price Index for All Urban Consumers: All Items",
    ylab= "Percent",
    source= "Data Source: Federal Reserve Bank of St.Louis (CPIAUCNS)",
    footnoteLeft  = "not seasonally adjusted",
    footnoteRight = paste("Extracted:", date()),
    lastObs = TRUE )
dev.off()

  yahoo <- TSconnect("getSymbol", dbname="yahoo") 
  x <- TSget("F", con=yahoo)
  #plot(x)

png(file=files[3],width=480, height=240, pointsize=12, bg = "white")
# mv Fordsmall.png Fordsmall.png.orig ; pngcrush Fordsmall.png.orig Fordsmall.png
#png(file="Ford.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(x[,"F.Close"], 
    Title = "Ford Motor Co. Closing Price", 
    ylab= "US$",
    source= "Data Source: yahoo symbol F",
    footnoteRight = paste("Extracted:", date()),
    lastObs = TRUE )
dev.off()

unlink(tmp, recursive = TRUE)
