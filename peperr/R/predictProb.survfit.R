predictProb.survfit <- function (object, response, x, times, train.data=NULL, ...){
   km.pred <- matrix(summary(object, times=times)$surv, nrow(x), length(times),byrow=TRUE)  
}