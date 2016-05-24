logis.noise <-
function(model,intensity=25) {
  lg<-nrow(model$model)
  int<-min(abs(model$resid))/intensity
  noise<-rnorm(lg,0,int)
  fit<-model$fitted.values
  values<-integer(lg)
  for (i in 1:lg) {
    if (fit[i]+noise[i]>1 | fit[i]+noise[i]<0) {
	values[i]<-fit[i]-noise[i]
    } else {
	values[i]<-fit[i]+noise[i]
    }
  }
  return(values)
}

