de <-
function(x, tree){
  if(length(x) != 1) stop("The length of x in function de must be 1.")
  y <- tree$node;  de <- NA
  if(sum(match(x, y), na.rm = TRUE) != 0) {
    temp <- 1:length(y)
    start <- match(x, y) + 1    
    end <- length(y)
    if(start <= length(y) & nchar(y[start]) > nchar(x)) {
      temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
      if(!is.na(temp1)) end <- temp1
      de <- y[start:end]
    }
  }
  de
}
