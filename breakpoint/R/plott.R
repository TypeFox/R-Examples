plott <-
function(data, x, L){
  ddx <- array(0, dim = c(L, 1))
  for(i in 1 : (length(x)- 1)){
    ddx[x[i] : (x[i + 1] - 1), ] <- (mean(data[x[i] : (x[i + 1] - 1), 1]))
  }
  return(ddx)
}
