"probhmm" <-
function (logalpha, logbeta, Pi, delta, cumprob) 
{
    n <- nrow(logalpha)
    prob <- rep(NA, n)
    for (i in 1:n) {
        if (i==1){
            pre <- delta
        }
        else {
            la <- logalpha[i-1,]
            pre <- exp(la - mean(la[la != -Inf])) %*% Pi
        }
        lb <- logbeta[i,]
        post <- exp(lb - mean(lb[lb != -Inf]))
        prob[i] <- (pre %*% diag(cumprob[i,]) %*% post)/(pre %*% post)
    }
    return(prob)
}
