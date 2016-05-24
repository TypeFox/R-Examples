myMatch <- function(p,x,tol=10e-6){

ind <- NA
  
i <- 1
while(i <= length(x)){
  if(abs(p-x[i])<tol){
    ind <- i
    i <- length(x)+1
  }
  else i <- i+1
}

return(ind)

}
