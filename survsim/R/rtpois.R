rtpois <-
function(N, lambda=NA)
{
   if (is.na(lambda) == TRUE) return(NA)
   qpois(runif(N, dpois(0, lambda), 1), lambda)
}
