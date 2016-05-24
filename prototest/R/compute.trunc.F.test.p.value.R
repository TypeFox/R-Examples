compute.trunc.F.test.p.value <-
function(F.stat, trunc.interval, df1, df2, verbose=FALSE){
  # construct the more extreme interval and intersect with the truncation interval
  # (which may be a union of distinct intervals)
  more.extreme.interval = Intervals(c(F.stat, Inf))
  if(verbose){
    print ('More extreme interval: ')
    print (more.extreme.interval)
    cat('\n')
    print ('Truncation interval: ')
    print (trunc.interval)
  }
  more.extreme.interval = interval_intersection(more.extreme.interval, trunc.interval)
  if(verbose){
    cat('\n')
    print ('After intersection with truncation interval: ')
    print (more.extreme.interval)
  }
  
  # compute the mass of the truncation interval (again, might be the union of disjoint intervals)
  trunc.interval.mass = sum(apply (trunc.interval, 1, function(row){
    pf(row[1], df1=df1, df2=df2, lower.tail=FALSE) - pf(row[2], df1=df1, df2=df2, lower.tail=FALSE)
  }))
  if (verbose){
    print (paste('Mass of truncation interval: ', trunc.interval.mass, sep=''))
  }
  
  # compute the mass of the more extreme interval intersected with the truncation interval
  more.extreme.interval.mass = sum(apply (more.extreme.interval, 1, function(row){
    pf(row[1], df1=df1, df2=df2, lower.tail=FALSE) - pf(row[2], df1=df1, df2=df2, lower.tail=FALSE)
  }))
  if (verbose){
    print (paste('Mass of more extreme interval: ', more.extreme.interval.mass, sep=''))
  }
  p.val = more.extreme.interval.mass/trunc.interval.mass
  if (verbose){
    print(paste('df1: ', df1, sep=''), quote=FALSE)
    print(paste('df2: ', df2, sep=''), quote=FALSE)
    print(paste('p-val: ', p.val, sep=''), quote=FALSE)
  }
  p.val
}
