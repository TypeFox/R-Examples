
stretch <- function(S,maxx){ 

  n = length(S$x) # a number of days in which survival function is estimated
  x = 0:maxx # days

  if (n==0) {
    y = array(1,length(x)) # if there is not any event, survival function is equal to function y=1
  } else {
    y = array(1,S$x[1]) # initial setting of a vector with survival estimates
    if (n>1) { 
      for(i in 1:(n-1)) {
        y=c(y,array(S$y[i],(S$x[i+1]-S$x[i]))) # assigns survival estimates to each day of the follow-up until the last step of the survival curve
      }
    }
    y=c(y,array(S$y[n],(maxx+1-S$x[n]))) # assigns survival estimates to each day between the last step of the survival curve and the maximum follow-up time
  }

  stretch  = list(x=x,y=y)
}
