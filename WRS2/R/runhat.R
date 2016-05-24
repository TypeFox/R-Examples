runhat<-function(x,y,pts=x,est=onestep,fr=1,nmin=1){
  #
  # running  interval smoother that can  be used  with any measure
  # of location or scale. By default, a modified one-step M-estimator is used.
  # This function computes an estimate of y for each x value stored in pts
  #
  # fr controls amount of smoothing
  rmd<-rep(NA,length(pts))
  for(i in 1:length(pts)){
    val<-y[near(x,pts[i],fr)]
    if(length(val)>=nmin)rmd[i]<-est(val)
  }
  rmd
}