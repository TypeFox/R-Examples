setClass("RFPHI", representation = representation(
    sectionName = "character",
    vectorMean = "numeric",
    harmonicMean = "numeric",
    strainRatio = "numeric",
    indexSymmetry = "numeric",
    sampleSize = "numeric",
    meanObjectArea = "numeric",
    rsAxes = "numeric",
    chiSquare = "data.frame"
  ))
  
setClass("FRY", representation = representation(
    sectionName = "character",
    vectorMean = "numeric",
    strainRatio = "numeric",
    sampleSize = "numeric",
    rsAxes = "numeric",
	meanObjectArea = "numeric",
    fryParams = "data.frame",
	voidScale = "numeric"
  ))
  