
#' @include 04ControlClass.R
NULL

#' defaultParameters
#' 
#' defaultParameters include nine experimental constructed features (techically, subclasses)
#' @export

defaultParameters <- initializeparameterclassobject(list("MissingValueShare", "MissingValueToClass", "LOFScore", "MahalanobisToClassCenter", "DistanceToNearest", "LenghtOfIQR", "NearestPointPurity", "ClassificationCertainty", "NeighborhoodDiversity", "ScatterCounter"))

