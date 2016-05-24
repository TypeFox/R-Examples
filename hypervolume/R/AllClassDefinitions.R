setClass("Hypervolume", slots=c(
    Name="character",
    Data="matrix",
    Dimensionality="numeric",
    Volume="numeric",
    PointDensity="numeric",
    Bandwidth="numeric",
    RepsPerPoint="numeric",
    DisjunctFactor="numeric",
    QuantileThresholdDesired="numeric",
    QuantileThresholdObtained="numeric",
    RandomUniformPointsThresholded="matrix",
    ProbabilityDensityAtRandomUniformPoints="numeric"
    ))

setClass("HypervolumeList", slots=c(
    HVList="list"
  ))


summary.Hypervolume <- function(object, ...)
{
  cat(sprintf("Hypervolume\n\tName: %s\n\tNr. of observations: %d\n\tDimensionality: %d\n\tVolume: %f\n\tBandwidth: %s\n\tDisjunct factor: %f\n\tQuantile desired: %f\n\tQuantile obtained: %f\n\tNr. of repetitions per point: %.0f\n\tNumber of random points: %.0f\n\tPoint density: %f\n", 
              object@Name, ifelse(all(is.nan(object@Data)), 0, nrow(object@Data)), object@Dimensionality, object@Volume, paste(format(object@Bandwidth,digits=2),collapse=' '), object@DisjunctFactor, object@QuantileThresholdDesired, object@QuantileThresholdObtained, object@RepsPerPoint, nrow(object@RandomUniformPointsThresholded), object@PointDensity))
  
}

summary.HypervolumeList <- function(object, ...)
{
  cat(sprintf("HypervolumeList with %d elements:\n\n", length(object@HVList)))
  
  if (length(object@HVList)>0)
  {
    for (i in 1:length(object@HVList))
    {
      whichhv <- object@HVList[[i]]
      
      summary.Hypervolume(whichhv)
    }
  }
}

setMethod("show","Hypervolume", function(object) {summary.Hypervolume(object)})
setMethod("show","HypervolumeList", function(object) {summary.HypervolumeList(object)})


get_volume.Hypervolume <- function(object)
{
  return(object@Volume)
}

get_volume.HypervolumeList <- function(object)
{
  sapply(object@HVList, get_volume.Hypervolume)
} 

hypervolume_join <- function(...)
{
  hvl <- list()
  
  for (a in list(...))
  {
    if (class(a) == "HypervolumeList")
    {
      hvl <- c(hvl, a@HVList)
    }
    else if (class(a) == "Hypervolume")
    {
      hvl <- c(hvl, a)
    }
  }
  
  return(new("HypervolumeList",HVList=hvl))
}

setGeneric("get_volume", function(object) {})
setMethod("get_volume","Hypervolume", function(object) {get_volume.Hypervolume(object)})
setMethod("get_volume","HypervolumeList", function(object) {get_volume.HypervolumeList(object)})
