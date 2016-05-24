ciLevel <-
function(x, n, p){
  if (x == 0){
    level <- c(0, p)
  } else if (x == n){
    level <- c(1 - p, 1)
  } else {
    level <- c((1 - p) / 2, 1 - (1 - p) / 2)
  }
  return(c(level, (1 - p) / 2))
}