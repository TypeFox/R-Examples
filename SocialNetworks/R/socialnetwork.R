

#' Generates social network based on xy spatial coordinates of individuals.
#' 
#' @param x a list of x coordinates for all the individuals  
#' @param y a list of y coordinates for all the individuals
#' @param ir a list of interaction radii for all the individuals
#'@description  This method uses the pairwise distances between each pair of individuals.
#'For a spatial point pattern X, the association of individual j on 
#'individual i, Aij,  is calculated using the distances between points X_{i} and X_{j} using a smooth interaction function
#'first introduced by Illian (2009).
#' If d represents the distance between points X_{i} and X_{j}, and the interaction radius for individual i is R_{i}, then the
#' association of j on i is calculated as:   ((1-(d/R_{i})^{2}))^{2} if d>0 and d<=R_{i}, and 0 otherwise.
#'This function has been described as a smooth interaction function because the value of the association calculated decreases smoothly as a function of
#'distance (See Figure (b)).  This is in contrast to associations calcuated using a binary function where the association =1 if d<=R and 0 otherwise.  Such a function is based
#'on the assumption that the association is constant (1) at all distances less than R_{i}, and 0 for distances greater than R_{i} (See Figure (a)).
#'
#'\figure{step.jpeg}
#'\figure{smooth.jpeg}
#'@references Illian, Moller, Waagepetersen, 2009. Hierarchical spatial point process analysis for a plant community of high biodiversity.Environ.
#'Ecol. Stat. vol 16, pp 389-405

#' @export
calculateassociations = function(x,y,ir){
  n = length(x)
  spnet = matrix(1,nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      d = sqrt( ((x[i] - x[j])^2)  +  ((y[i]-y[j])^2) )
      #spnet[i,j] =  ifelse( d>0&&d<ir[i], (1-(d/ir[i])^2)^2,0)
      spnet[i,j] =  ifelse( ((d>0)&&(d<=ir[i])), (1-(d/ir[i])^2)^2,0)
    }
  }
  
  return(spnet)
}


