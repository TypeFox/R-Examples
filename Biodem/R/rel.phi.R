rel.phi <- function(x,R,method="A"){
  metodi <- c("A","B")
  method <- pmatch(method, metodi)
  if (is.na(method))
    stop("not valid method")
  if (method==1){
    x1 <- x/4
    x1
  }
  else{
    x1 <- x/4+(3*R*((x-R))/(16*(1-R)))
    x1
  }
}
