F.double.obs.prob <- function( df, observer = "both" ){
#
#   Compute the probability of detection from a double observer system.
#   No external covariates allowed here.
#
#   Inputs:
#   df = data frame containing the components $obsby.1, $obsby.2.
#       These components are TRUE/FALSE (logical) vectors indicating whether
#       observer 1 (obsby.1) or observer 2 (obsby.2) spotted the target.
#   observer = indicates whether observer 1 or observer 2 or both were full-time observers.
#       If, for example, observer 2 was a data recorder and part-time observer, or if observer 2
#       was the pilot, set "observer = 1".  This dictates which set of observations form the denominator
#       of the double observer system.  For example, if "observer = 1", observations by observer 1 that were not seen
#       by observer 2 are ignored. The estimate in this case uses targets seen by both observers and
#       those seen by observer 2 but not observer 1. If observer = "both", the computation goes both directions.

if( !(all(c("obsby.1", "obsby.2") %in% names(df)))){
    stop("Variables 'obsby.1' and 'obsby.2' not found in input data frame.")
}

obs1 <- as.logical( df$obsby.1 )
obs2 <- as.logical( df$obsby.2 )
obs.both  <-  obs1 & obs2

if( is.character( observer ) ){
    if(observer == "both"){
        obs.tot  <-  obs1 | obs2  # this should be all 1's, assuming no extra lines are input (like "unknown")
        p1 <- sum( obs.both ) / sum( obs2 )
        p2 <- sum( obs.both ) / sum( obs1 )
        # Assume observers are independent here.  This was checked via simulation. This estimator is close to unbiased when observers are actually independent.
        p <- p1 + p2 - p1*p2
    } else {
        stop("Inappropriate 'observer' parameter")
    }
} else if( observer == 1 ){
    p <- sum( obs.both ) / sum( obs2 )
} else if( observer == 2 ){
    p <- sum( obs.both ) / sum( obs1 )
} else {
    stop("Inappropriate 'observer' parameter")
}

p

}        
        
