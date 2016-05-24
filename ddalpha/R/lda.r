lda.train <- function(data){
  new.frm <- data.frame(data)
  z <- lda(formula=as.formula(paste("X", ncol(data), " ~ .", sep="")), data = new.frm, tol = sqrt(.Machine$double.eps))#, prior = rep(1, dimension)/dimension)
  return (z)
}

lda.classify <- function(objects, lda){
  z <- predict(lda, data.frame(objects))$class
  return (z)
}
