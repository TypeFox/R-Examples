# define pmf (vectorized), 
# checking to be sure input is in { 0, 1, 2, 3, 4 }
f <- function(x) {
    sapply(x, function(x) {
        if ( ! ( x %in% 0:4 ) ) { return(0) }
        return( factorial(4) / ( 16 * factorial(x) * factorial(4-x) ) )
    })
}
f(0:6)   
probplot1 <- xyplot(f(0:4)~0:4, xlab="x", ylab="probability")
probplot2 <- xyplot(f(0:4)~0:4, xlab="x", ylab="probability", 
                    type="h")
probplot3 <- xyplot(f(0:4)~0:4, xlab="x", ylab="probability", 
                    type=c("l","p"))
