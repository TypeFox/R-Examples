#' Extrapolate Horvitz-Thompson abundance estimates to entire surveyed region
#'
#' @param Nhat.by.sample dataframe of abundance by sample
#' @param samples samples table
#' @param width transect width
#' @param point if TRUE point count otherwise line transect
#' @return Revised Nhat.by.sample dataframe containing estimates extrapolated
#'   to survey region
#' @note Internal function called by \code{\link{dht}} and related functions.
#' @author Jeff Laake
#' @keywords utility
survey.region.dht <- function(Nhat.by.sample, samples, width, point){
  #  Compute effort in each region and the area in the covered region
  Effort.by.region <- by(samples$Effort,samples$Region.Label,sum)
  if(point){
    CoveredArea <- pi*as.vector(Effort.by.region)*width^2
  }else{
    CoveredArea <- 2*as.vector(Effort.by.region)*width
  }

  # Scale up abundance in covered region to the survey region
  # unless no areas given
  Nhat.by.sample <- merge(Nhat.by.sample,
                          data.frame(Region.Label=names(Effort.by.region),
                                     CoveredArea=CoveredArea,
                                     Effort=as.vector(Effort.by.region)),
                          by.x="Region.Label",
                          by.y="Region.Label",
                          all.x=TRUE)
  Nhat.by.sample$Nhat <- Nhat.by.sample$Nhat*Nhat.by.sample$Area/
                          Nhat.by.sample$CoveredArea
  return(Nhat.by.sample)
}
