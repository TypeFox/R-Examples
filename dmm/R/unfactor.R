unfactor <-
function(x){
# unfactor() - convert vector x from factor to numeric
#              non numeric values will coerce to NA
  if(!is.factor(x)){
    stop("unfactor: argument x must be a factor:\n")
  }
  y <- as.numeric(as.character(x))
  return(y)
}
