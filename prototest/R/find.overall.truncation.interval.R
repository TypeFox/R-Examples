##### finds the overall interval (intersection of component ones) associated the Q, R and S vectors
##### input:
#####   - Q, R, S = Vectors of the same length. They contain replications of the q, r and s parameters respectivelys
find.overall.truncation.interval <-
function(Q, R, S, verbose=FALSE){
  num.intervals = length (Q)
  if (verbose){print (paste('Number of intervals to compute: ', num.intervals), quote=FALSE)}
  
  #overall.interval = Intervals (c(0, Inf))
  #overall.interval =lapply(1:num.intervals, function(i){    
  overall.interval = do.call(interval_intersection, lapply(1:num.intervals, function(i){    
    if(verbose){
      print ('------------------------------------------')
      print (paste('Interval ', i, sep=''), quote=FALSE)
      print ('', quote=FALSE)
    }
    interval = find.single.truncation.interval(Q[i], R[i], S[i], verbose=verbose)
    if (verbose){print (interval)}
    interval
  }))
  
  overall.interval
}
