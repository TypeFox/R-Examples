paintPoints <- function(aDetectedPlot, aBasePlot, size = 8, color = "blue") {
  
  blankImage <- channel(Image(dim = dim(aDetectedPlot)),  "rgb")
  theCoordinates <- computeFeatures.moment(aDetectedPlot)[, 1:2]
  kern = makeBrush(size = size, shape = "disc", step = FALSE)
  median <- median(1:size) - 1 - 0.5
  
  for(a in 1:nrow(theCoordinates)) {
    corX <- (theCoordinates[a, 1] - median)
    corY <- (theCoordinates[a, 2] - median)
    blankImage[corX:(corX + size - 1 + 0.5), corY:(corY + size - 1), 1] <- kern
  }
    
  blankImage <- channel(blankImage, "grey")
  
  theNew <- paintObjects(blankImage, 
                         channel(aBasePlot, "rgb"), 
                         col = rgb(t(col2rgb(color)), maxColorValue = 255), 
                         opac = 0.5,
                         thick = TRUE) 
  
  return(theNew)
}

unitTrim <- function(unit) {
  as.numeric(sub("native", "", as.character(unit)))
}