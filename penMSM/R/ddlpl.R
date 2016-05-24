ddlpl <- function(b, X, Ri, Ci){
    n <- nrow(X)
    ef <- exp(as.numeric(X%*%b))
    e2f <- exp(2*as.numeric(X%*%b))
    dummy <- dummy2 <- secondpart <- secondpart2 <- rep(0, n)
    for (i in 1:n) {
        dummy[i] <- sum(ef[Ri[[i]]])
        dummy2[i] <- dummy[i]*dummy[i]}
    dummy[which(dummy==0)] <- 1e-05
    dummy2[which(dummy2==0)] <- 1e-05
    for (i in 1:n) {
        secondpart[i] <- sum(1/dummy[Ci[[i]]])
        secondpart2[i] <- sum(1/dummy2[Ci[[i]]])}
    secondgradients <- e2f*secondpart2 - ef*secondpart
    return(secondgradients)}
