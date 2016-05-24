pPool <- function(p){
    n <- length(p)
    k <- 1 + floor(n / 2)
    pp <- rev(sort(p))
    
    j <- 2:k
    red <- c(1, pp[2 * j - 2])
    
    return(red)
}
