

#' Makes a social network based on the xy spatial coordinates provided.
#' 
#' @param arg1 x coordinates for individuals
#' @param arg2 y coordinates for individuals
#' @param arg3 interaction radii for each individual (they can all be equal)
#' @param numpts number of Monte Carlo simulations

#' #'@description  This method uses the area of the intersection between the territorial zones
#' or home bases between each pair of individuals.  For a spatial point pattern X, the association of individual j on 
#'individual i, Aij,  is calculated as the percentage of the overlap of the discs centered at 
#' points X_{i} and X_{j} of the total area of the territorial area for individual i.  The radius for each disc is the inputted interaction radius.
#' The interaction radius for a given population can be identical for each individual, or
#' different.  The interaction radius represents the area within which an individual extracts
#' nutrients or exerts its influence, or communicates an action.
#' 
#' 
#' The associations calcuated using this method can be asymmetric.  In this case, the
#' interaction radii for two given individuals would be different, implying that the 
#' proportion of the overlap between the zones for the individuals is different for each 
#' individual. As as example, Figure 1 illustrates the effect of different interaction radii
#' per individual.  Individual i is represented by the filled square and individual j is represented
#' by the filled circle. The percentage of the overlap between the two territorial zones in the total
#' area of territorial zone i is larger than that in territorial j, suggesting that the effect of 
#' individual j on i is greater than that of i on j. 
#' 
#' The calculations are done based on a Monte Carlo method.
#'
#'\figure{asymmetry.jpeg}

#' @export



calculate.areas <- function(arg1, arg2, arg3, numpts) .Call("rcpp_calculate_areas", arg1, arg2, arg3, runif(numpts, 0, 1), runif(numpts, 0, 1), PACKAGE="SocialNetworks")