
predict.clogitboost <- function(object, x, strata, ...){
  x <- as.matrix(x)
  obj <- fitclogitboostpredict(object, x, strata)

  prediction <- matrix(0, ncol = 1, nrow = length(obj$prob))
  for (i in seq(1, length(obj$prob))) {
    m <- max(obj$prob[which(strata == strata[i])])
    if (obj$prob[i] == m) {
      prediction[i] <- 1
    } else {
      prediction[i] <- 0
    }
  }
  
  return(list(prob = obj$prob, utility = obj$utility, prediction = prediction))
}      
