xinner_kernel <- function(x, y, kernel, kparam){

  K.train <- Kmat(x = x, 
                  y = x, 
                  kernel = kernel, 
                  kparam = kparam)

  xinner <- K.train + 1.0

  return(xinner)

}

xinner_linear <- function(x){
  xinner <- x %*% t(x)
  return(xinner)
}
