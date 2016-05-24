qda.train <- function(data){
  new.frm <- data.frame(data)
  z <- qda(formula=as.formula(paste("X", ncol(data), " ~ .", sep="")), data = new.frm)#, prior = rep(1, dimension)/dimension)
  return (z)
}

qda.classify <- function(objects, qda){
  z <- predict(qda, data.frame(objects))$class
  return (z)
}
