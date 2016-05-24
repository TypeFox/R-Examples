eql <- function(x,y) {
  if (length(x)>length(y)) y <- rep(y,,length(x)) else x <- rep(x,,length(y))
  ifelse(is.na(x),is.na(y),ifelse(is.na(y),FALSE,x==y))
}
