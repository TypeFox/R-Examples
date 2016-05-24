seq2 <-
function(x,int=999) {
  x2 <- na.omit(x)
  result <- seq(min(x2),max(x2),abs(max(x2)-min(x2))/int)
  return(result)
}

