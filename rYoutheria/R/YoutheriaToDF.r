#' Convert data returned from API to a dataframe
#' 
#' Takes the list returned by \code{fromJSON(getURL(url))}, where
#' url is the ValueByType controller in the YouTheria API, and returns
#' a data.frame. This code is faster than using plyr.
#' 
#' @param x a list of trait data as returned by \code{fromJSON(getURL(url))}.
#' 
#' @return A \code{data.frame} of melted trait data        

YoutheriaToDF <- function(x){ 
  
  MeasurementTypeID <- vapply(x, function(x) as.numeric(x$MeasurementTypeID),0)
  MeasurementSetID <- vapply(x, function(x) as.numeric(x$MeasurementSetID),0)
  StudyUnitId <- vapply(x, function(x) as.numeric(x$StudyUnitId),0)
  Genus <- vapply(x, function(x) as.character(x$Genus),"")
  Species <- vapply(x, function(x) as.character(x$Species),"")
  SubSpecies <- vapply(x, function(x) as.character(x$SubSpecies),"")
  MSW93Binomial <- vapply(x, function(x) as.character(x$MSW93Binomial),"")
  MSW05Binomial <- vapply(x, function(x) as.character(x$MSW05Binomial),"")
  AuthorityText <- vapply(x, function(x) as.character(x$AuthorityText),"")
  ValueType <- vapply(x, function(x) as.character(x$ValueType),"")
  MValue <- vapply(x, function(x) as.character(x$MValue),"")
  
  r <- data.frame(MeasurementTypeID,MeasurementSetID,StudyUnitId,Genus,Species,SubSpecies,
                  MSW93Binomial,MSW05Binomial,AuthorityText,ValueType,MValue)
  return(r)
}