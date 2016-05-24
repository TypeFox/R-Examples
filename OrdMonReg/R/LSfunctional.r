LSfunctional <- function(f1, g1, w1, f2, g2, w2){

    res <- sum(w1 * (f1 - g1) ^ 2) + sum(w2 * (f2 - g2) ^ 2)
    return(res)
}
