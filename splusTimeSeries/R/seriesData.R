"seriesDataNew" <- 
function()
matrix(nrow=0, ncol=0)

"seriesData" <- 
function(object)
{
  ## return the data inside an ordered data object
  object@data
}

"seriesData<-" <- 
function(object, value)
{
  ## replace the data inside an ordered data object
  value <- asSeriesData(value)
  if(length(object@positions) != numRows(value))
    stop("Positions and data lengths do not agree")
  object@data <- value
  object
}

"seriesDataValid" <- 
function(object)
is.rectangular(object)

"asSeriesData" <- 
function(object)
as.rectangular(object)



