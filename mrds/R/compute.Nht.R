#' Horvitz-Thompson estimates 1/p_i or s_i/p_i
#'
#' Compute individual components of Horvitz-Thompson abundance estimate in
#' covered region for a particular subset of the data depending on value of
#' group = TRUE (do group abundance); FALSE(do individual abundance)
#'
#' @param pdot vector of estimated detection probabilities
#' @param group if TRUE (do group abundance); FALSE(do individual abundance)
#' @param size vector of group size values for clustered populations
#' @return vector of H-T components for abundance estimate
#' @note Internal function called by \code{\link{covered.region.dht}}
#' @author Jeff Laake
#' @keywords utility
compute.Nht <- function(pdot,group=TRUE,size=NULL){
  if(group){
     return(1/pdot)
  }else{
    if(!is.null(size)){
      if(!any(is.na(size))){
        return(size/pdot)
      }else{
        stop("One or more missing group sizes in observations")
      }
    }else{
      stop("Missing group sizes in observations")
    }
  }
}
