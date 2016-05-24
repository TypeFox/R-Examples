################ stored matrices from http://www2.research.att.com/~njas/hadamard/ 
################ filling the gaps up to 256, 260 is next gap

hadamard.table <- function(n){
  if(n %in% c(52,92,100,116,156,172,188,236,244)){
    data(hadamard.table,envir = environment())
    H<-get(paste("had",n,sep=""))
  } else {
    warning(paste("no hadamard matrix stored for this size:",n))
    H <- NULL
  }
  H
}
