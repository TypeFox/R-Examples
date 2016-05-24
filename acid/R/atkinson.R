atkinson <-
function (x, epsilon = 1) {
    x <- as.numeric(x)
    if (is.null(epsilon)) epsilon <- 1
    if (epsilon == 1){
        A <- 1 - (exp(mean(log(x)))/mean(x))
    }else{
        x <- (x/mean(x))^(1 - epsilon)
        A <- 1 - mean(x)^(1/(1 - epsilon))
    }
    return(A)
}
