sample.tte <-
function(time,surv) {
  m=length(time)
  n=ncol(surv)
  if(m!=nrow(surv)) stop("Mismatch time and surv!")
  times=rep(time[m],n)
  event=rep(0,n)
  for(i in 1:n) {
    u=runif(1)
    for(j in 1:m) {
      if(u>surv[j,i]) { times[i]=time[j]; event[i]=1 ; break } 
    }
  }
  return(list(times=times,event=event))
}
