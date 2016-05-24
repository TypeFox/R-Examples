myUnique <- function(x,tol=10e-6){

n <- length(x)

#print(x)

i <- 1
while(i <= length(x)){
  ind <- abs(x[i]-x) > tol
  #print(ind)
  x <- c(x[i],x[ind])
  #print(x)
  i <- i+1
}

return(x)

}
