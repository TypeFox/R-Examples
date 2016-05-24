rrice <- function(n, vee, sigma){

    if ((!is.numeric(n)) || (length(n)>1)){ 
        stop("bad input for argument 'n'")
    }
    theta <- 1
    X <- rnorm(n, mean = vee * cos(theta), sd = sigma)
    Y <- rnorm(n, mean = vee * sin(theta), sd = sigma)
    sqrt(X^2 + Y^2)*sqrt(2)
}
