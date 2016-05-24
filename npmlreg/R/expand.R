"expand" <-
function(x,k){
  xx <- x
  for ( i in 2:k) xx <- rbind(x,xx)
  xx}

