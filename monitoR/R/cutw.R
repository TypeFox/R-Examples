# For cutting a wave object down to specified time limits
# Modified: 2015 August 1

cutw <-
function(
   wave,       # Wave object to be cut
   from=NULL,  # Initial time to trim to (s)
   to=NULL     # Final time to trim to (s)
) {

  .Deprecated(cutWave, package = "monitoR", msg = "The cutw() function in the monitoR package has been renamed cutWave() \nto avoid conflict with the cutw() function in the seewave package.\nIf you want the monitoR version, please switch to cutWave().\nIf you want the seewave version, load seewave after monitoR or \"use seewave::cutw()\".\nThe monitoR version of cutw() will be removed in the next update. " )

   if(is.null(from)) from <- 0
   if(is.null(to)) to <- length(wave@left)/wave@samp.rate
   i.from <- floor(from*wave@samp.rate)
   i.to <- ceiling(to*wave@samp.rate)
   wave@left <- wave@left[i.from:i.to]
   if(wave@stereo) wave@right <- wave@right[i.from:i.to]

   return(wave)
}
