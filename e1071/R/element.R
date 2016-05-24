element <- function(x, i) {

  if(!is.array(x))
    stop("x is not an array")
  
  ni <- length(i)
  dx <- dim(x)
  
  if(length(i)!=length(dx))
    stop("Wrong number of subscripts")

  if(ni==1){
    return(x[i])
  }
  else{
    m1 <- c(i[1], i[2:ni]-1)
    m2 <- c(1,cumprod(dx)[1:(ni-1)])
    return(x[sum(m1*m2)])
  }
}

