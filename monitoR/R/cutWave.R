# For cutting a wave object down to specified time limits
# Modified: 2015 August 1

cutWave <-
function(
   wave,       # Wave object to be cut
   from=NULL,  # Initial time to trim to (s)
   to=NULL     # Final time to trim to (s)
) {

   if(is.null(from)) from <- 0
   if(is.null(to)) to <- length(wave@left)/wave@samp.rate
   i.from <- floor(from*wave@samp.rate)
   i.to <- ceiling(to*wave@samp.rate)
   wave@left <- wave@left[i.from:i.to]
   if(wave@stereo) wave@right <- wave@right[i.from:i.to]

   return(wave)
}
