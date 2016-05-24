"wvec" <- 
function(x,y = y) {
   w <- NULL
   w[y == 1] <- 1-x
   w[y ==-1] <- x
   w
}
