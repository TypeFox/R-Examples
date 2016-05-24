lpl <- function(beta, X, risksetlist, event){
    n <- length(event)
    risk <- rep(0, n)
    f <- as.numeric(X%*%beta)
    ef <- exp(f)
    for (i in which(event == 1)){## 1:n) {
        risk[i] <- sum(ef[risksetlist[[i]]])
    }
    logpartiallikelihood <- sum(event * (f - log(risk)))
    return(logpartiallikelihood)}