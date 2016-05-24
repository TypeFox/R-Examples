dlpl <- function(event, b, X, Ri, Ci){
    n <- length(event)
    ef <- exp(as.numeric(X%*%b))
    dummy <- secondpart <- rep(0, n)
    for (i in 1:n) {
        dummy[i] <- sum(ef[Ri[[i]]])
    }
    dummy[which(dummy==0)] <- 1e-05
    for (i in 1:n) {
        secondpart[i] <- sum(1/dummy[Ci[[i]]])
    }
    gradients <- event - (ef*secondpart)
    return(gradients)}
