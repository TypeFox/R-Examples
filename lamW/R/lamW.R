lambertW0 <- function(x){
  LAM <- double(length(x))
  LAM <- lambertW0_C(x)
  return(LAM)
}

lambertWm1 <- function(x){
  LAM <- double(length(x))
  LAM <- lambertWm1_C(x)
  return(LAM)
}