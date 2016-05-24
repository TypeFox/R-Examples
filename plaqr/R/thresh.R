
threshold <- function(fit, t, newdata=NULL, ...)
{
  pred <- predict(fit, ...)
  train.class <- 0*(pred < t) + 1*(pred >= t)
  true.class <- 0*(fit$y < t) + 1*(fit$y >= t)

  difference <- train.class - true.class
  train.error <- mean(difference != 0)

  true.low <- sum(true.class==0)
  true.high <- sum(true.class==1)
  false.low <- sum(difference == -1)
  false.high <- sum(difference == 1)

  pred.class <- NULL

  if(!is.null(newdata)){
    pred <- predict(fit, newdata=newdata, ...)
    pred.class <- 0*(pred < t) + 1*(pred >= t)
  }

  thresh <- list(pred.class=pred.class, t=t, train.class=train.class,
                 train.error=train.error, true.class=true.class,
                 true.high=true.high, true.low=true.low,
                 false.high=false.high, false.low=false.low,
                 call=fit$call, formula=fit$formula)
  class(thresh) <- c("thresh", "list")
  return(thresh)
}



print.thresh <- function(x,...)
{
  cat("Regression Call:\n")
  print(x$call)
  cat("\nThreshold: ",x$t,"\n")
  cat("Prediction: ",!is.null(x$pred.class),"\n")
  cat("True High: ",x$true.high,"    False Low:  ",x$false.low,"\n")
  cat("True Low:  ",x$true.low, "    False High: ",x$false.high,"\n")
}