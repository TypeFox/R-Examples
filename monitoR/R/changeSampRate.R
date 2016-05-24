# Modified: 2015 JULY 21

changeSampRate <- function(
   wchange,                # Wave object, sampling rate will be changed 
   wkeep=NULL,             # Wave object that has the sampling rate needed 
   sr.new=wkeep@samp.rate, # Or just specify the new sampling rate (Hz)
   dither=FALSE,           # Dither resampled wave object with gaussian noise
   dith.noise=32           # Std Dev of guassian noise function; 7 bits error (2^7=+/-2sd) adds ~ 6 dB noise
   ) {

   f <- wchange@samp.rate/sr.new
   if(dither) {
      dither <- rnorm(floor(length(wchange@left)/f), 0,dith.noise)
      wchange@left <- round(spline(x=1:length(wchange@left), y=wchange@left, xout=1:floor(length(wchange@left)/f)*f)$y+dither)
      if(wchange@stereo)
         wchange@right <- round(spline(x=1:length(wchange@right), y=wchange@right, xout=1:floor(length(wchange@right)/f)*f)$y+dither) 
   } else {  
      wchange@left <- round(spline(x=1:length(wchange@left), y=wchange@left, xout=1:floor(length(wchange@left)/f)*f)$y)
      if(wchange@stereo)
         wchange@right <- round(spline(x=1:length(wchange@right), y=wchange@right, xout=1:floor(length(wchange@right)/f)*f)$y)
   } 
   wchange@samp.rate <- sr.new

   return(wchange)
}
