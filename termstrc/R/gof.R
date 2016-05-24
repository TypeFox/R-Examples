
## Root mean squared error
rmse <- function (actual,estimated) {
  e <- actual - estimated
  sqrt(mean(e^2))
}

## Average absolute error
aabse <-function (actual,estimated){
  e <- actual - estimated	
  mean(abs(e))
}   							



