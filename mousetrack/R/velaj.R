## [veloc accel jerk] = velaj(t,x,y,winsz)
## time in milliseconds
## computes smoothed velocity, acceleration, jerk, from trajectory

.packageName <- 'mousetrack'

velaj <- function(t, x, y, winsz){

  tdiff = t[winsz:length(t)] - t[1:(length(t) - (winsz-1))]

  ##  subtract first pair
  ydiff = y[1:(length(y)-(winsz-1))] - y[2:(length(y)-(winsz-2))];
  xdiff = x[1:(length(x)-(winsz-1))] - x[2:(length(x)-(winsz-2))];
  dist = sqrt(ydiff^2 + xdiff^2);

  ##  add next pairs w/i window
  for (indx in 2:(winsz-1)){
    
    ydiff = y[indx: ( length(y)-(winsz-indx) ) ] -
      y[(indx + 1): (length(y) - (winsz-indx-1)) ]
    
    xdiff = x[indx: (length(x) - (winsz-indx)) ] -
      x[(indx + 1): (length(x) - (winsz-indx-1)) ]
    
    dist = dist + sqrt(ydiff^2 + xdiff^2) 
  }
  
  veloc = (dist/tdiff)*1000
  accel = veloc[3:length(veloc)] - veloc[1: (length(veloc)-2)]
  jerk =  accel[3:length(accel)] - accel[1: (length(accel)-2)]
  
# dist=sqrt(ydiff.^2+xdiff.^2);
# vel=(dist./tdiff)*1000;
# acc=vel(2:length(vel))-vel(1:length(vel)-1);
# jerk=acc(2:length(acc))-acc(1:length(acc)-1);

  return( list( veloc = veloc, accel = accel, jerk = jerk) )
       
}
