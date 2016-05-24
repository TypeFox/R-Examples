ftilde <- function(loc, yy, fhat, dy, h, type){

    # f0tilde
    if (identical(type, 0)){res <- sum(dnorm(yy, loc, h) * fhat) * dy}
    
    # f1tilde
    if (identical(type, 1)){res <- sum(-(loc - yy) * dnorm(yy, loc, h) * fhat) * dy / (h ^ 2)}
    
    # f2tilde
    if (identical(type, 2)){res	<- sum((((loc - yy) / h) ^ 2 - 1) * dnorm(yy, loc, h) * fhat) * dy / (h ^ 2)}

    return(res)
}
