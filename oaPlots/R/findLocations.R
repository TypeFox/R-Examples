#' Returns a Vector of x Locations
#' @param n number of observations for a given value
#' @param space space between points
#' @param center center plotting value 
#' @return numeric vector of location values
#' @author Jason Waddell
#' @export
findLocations <- function(n, space, center){
  
  locations <- rep(0, n)
  odd <- n %% 2
  
  if(odd){
    locations[1] <- center
    if(n > 1)
      for(i in 2:n)
        locations[i] <- ifelse(i %% 2 == 0, 
            center-floor(i/2)*space,
            center+floor(i/2)*space)
  }
  
  if(!odd){
    for(i in 1:n)
      locations[i] <- ifelse(i %% 2 == 0, 
          center-ceiling(i/2)*space + space/2,
          center+ceiling(i/2)*space - space/2)
  }
  return(locations)
}

