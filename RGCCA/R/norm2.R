norm2 <-
function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}
