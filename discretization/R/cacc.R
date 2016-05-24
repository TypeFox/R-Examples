cacc <-
function(tb){
        n <- sum(tb)
        nr <- dim(tb)[1]+1
         logn <- 0
         if(nr>1)logn <- log(nr)
         yp <- chiSq(tb)
         val=sqrt(yp/(yp+n*logn))
         return(val)
}
