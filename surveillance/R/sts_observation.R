################################################################################
#' Function for creating a sts-object with a given observation date
################################################################################
# Parameters 
###
#' @param sts sts-object we want to set at a previous state. Needs to include a reporting triangle.
#' @param dateObservation Date for which we want the state. Needs to be in the reporting triangle dates.
#' @param cut Boolean indicating wether to have 0 counts after the observation date or to simply cut the sts-object
#' @examples
#' data("salmAllOnset")
#' salmAllOnsety2013m01d20 <- sts_observation(salmAllOnset,
#' dateObservation="2014-01-20",cut=FALSE)
#' plot(salmAllOnset)
#' lines(salmAllOnsety2013m01d20@@observed,t="h",col="red")
#' @export
sts_observation <- function(sts,dateObservation,cut=TRUE){
  # The sts object we shall return
  stsSub <- sts
  
  # Index of the observation date
  line1 <- which(epoch(sts)==dateObservation)
  
  # Maximal delay
  D <- dim(stsSub@control$reportingTriangle$n)[2]-1
  
  # Number of dates
  theEnd <- dim(stsSub@control$reportingTriangle$n)[1]
  
  # Nothing observed after the observation date (I am a genius)
  stsSub@control$reportingTriangle$n[(line1+1):theEnd,] <- NA
  stsSub@observed[(line1+1):theEnd] <- 0
  
  # Not everything observed before the observation date
  
  for (i in 1:D){
    stsSub@control$reportingTriangle$n[line1+1-i,(i+1):(D+1)] <- NA
    stsSub@observed[line1+1-i] <- sum(stsSub@control$reportingTriangle$n[line1+1-i,],na.rm=T)
  }
  stsSub@control$reportingTriangle$n <- stsSub@control$reportingTriangle$n[1:line1,]
  # Return the new sts object
  if (cut){return(stsSub[1:line1])}
  else{return(stsSub)}
  
  
}