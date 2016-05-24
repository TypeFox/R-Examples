globalVariables("places", package="birdring", add=TRUE)

place2name <-
function (x) 
{
  data(places, envir = environment())
  placesx <- places[match(as.character(x), places$code), c(1, 2, 4), drop=TRUE]
  placesx$country <- as.factor(placesx$country)
  placesx$region <- as.factor(placesx$region)
  rownames(placesx) <- seq(1:nrow(placesx))
  names(placesx) <- c('country.name', 'region.name', 'current')
  
  if( any(is.na(placesx$country)) ){ # check that all codes were assigned
    nmiss <- x[is.na(placesx$country)]
    warning(paste(length(nmiss), 'place codes not recognised:', 
                  paste(unique(nmiss), collapse=',')), call.=FALSE)
  }
  
  return(placesx)
}