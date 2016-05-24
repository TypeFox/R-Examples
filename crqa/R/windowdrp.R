## GNU License, written by Moreno I. Coco
## original Matlab code by Rick Dale

## calculate recurrence over different sized windows
## arguments: step = interval considered on the serie;
##            window_size = the size of the window wherin crqa is runned. 
##            lag_profile_width = lags within the window

# step = 10;
# window_size = 100;
# lag_profile_width = 20;

# source("drpd_fromts.R")

.packageName <- 'crqa'

windowdrp <- function(x, y, step, windowsize, lagwidth, datatype, radius){
  
x = as.vector(as.matrix(x));   y = as.vector(as.matrix(y));
points = seq(1, (length(x) - (windowsize)-1), step)

drpd = vector() ## a vector to store all the points

for (i in points){
  
  xwin = x[i:(i+windowsize)];
  ywin = y[i:(i+windowsize)];
  
  drpd = c(drpd,
      mean(drpdfromts(xwin, ywin, lagwidth, datatype, radius)$profile))
      ## can also take the max here

}

maxrec = max(drpd);
maxlag = which(drpd == maxrec);

return( list(profile = drpd, maxrec = maxrec, maxlag = maxlag) )

}


