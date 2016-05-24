brokenstick <-
function(a){
  q <- matrix(1,a,a)
  q[which(upper.tri(q,diag=T)==F)] <- 0 
  q <- q%*%(1/1:a)/a
  return(q)
}
