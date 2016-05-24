rel.cond <- function(x,R,method="A"){
  metodi <- c("A","B")
  method <- pmatch(method, metodi)
  if (is.na(method))
    stop("not valid method")
  if (method==1){
    x1 <- (x-R)/(4-R)
    x1
  }
  else{
    x1 <- (x-R)/(4*(1-R))
    x1
  }
}
