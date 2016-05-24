kmPredictWrapper <- function(newdata, km.object) 
  predict(object = km.object, newdata = newdata, type = "UK", 
          se.compute = FALSE, checkNames = FALSE)$mean
