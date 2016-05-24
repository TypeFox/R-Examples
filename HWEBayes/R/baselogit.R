baselogit <-
function(probs){
    k <- length(probs)
    baselogit <- log(probs[-k]) - log(probs[k])
    list(baselogit=baselogit)
}

