hramp.prime <-
function(x, lower, eps)
{
    if(length(lower) > 1){
        mapply(hramp.prime, x, lower, eps)
    } else{
        if (lower==-Inf){
            return(1)
        } else{
            dhr <- (x-lower)/eps
            dhr[x>(lower+eps)] <- 1
            dhr[x<lower] <- 0
            return(dhr)
        }
    }
}
