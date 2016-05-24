#' Suggested optimal value for N
#' 
#' Based on previous similutation work, the function suggests an optimal
#' value for the N parameter based on national prevalence, the total number
#' of observations and the number of clusters. See Larmarange et al. 2011
#' for more details.
#' 
#' @param object object of class \code{\link[=prevR-class]{prevR}}.
#' 
#' @return an integer.
#' 
#' @references
#' Larmarange Joseph, Vallo Roselyne, Yaro Seydou, Msellati Philippe and Meda Nicolas (2011) 
#' "Methods for mapping regional trends of HIV prevalence from Demographic and Health Surveys (DHS)", 
#' \emph{Cybergeo: European Journal of Geography}, no 558, \url{http://cybergeo.revues.org/24606}, 
#' DOI: 10.4000/cybergeo.24606.
#' 
#' @examples 
#' Noptim(fdhs)
#' 
#' @export
#' @keywords stat

Noptim <- function(object) {
  clusters               = slot(object,"clusters")
  clustersNumber         = nrow(clusters)
  ObservationNumber      = sum(clusters$n)
  PositiveCases          = sum(clusters$pos)
  isWeightedData         = !any(is.na(match(c("wn","wpos"),names(clusters))))
  nationalPrev           = 100*sum(clusters$pos,na.rm=T)/sum(clusters$n,na.rm=T)
  if(isWeightedData){
    nationalPrev = 100*sum(clusters$wpos,na.rm=T)/sum(clusters$wn,na.rm=T)
  }
  if (nationalPrev>50) {
    nationalPrev <- (100 - nationalPrev)
  } 
  round(14.172*ObservationNumber^0.419*nationalPrev^-0.361*clustersNumber^0.037-91.011)
}