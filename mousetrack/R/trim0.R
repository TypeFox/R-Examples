## trim the trajectory

.packageName <- 'mousetrack'

trim0 <- function(x, y, thresh){

  o1 = x[1]
  o2 = y[1]

  dists = sqrt((x-o1)^2 + (y-o2)^2)
  
  latindx = which( dists > thresh)[1] # latency index

  trimmed1 = x[latindx:length(x)];
  trimmed2 = y[latindx:length(y)];

return (list(latindx = latindx, xtrim = trimmed1, ytrim = trimmed2) )

}

