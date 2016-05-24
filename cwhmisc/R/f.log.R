f.log <- function(x){
    i <- (x>0)
    aus <- x[i]
    const <- median(aus)/((median(aus)/quantile(aus, probs = 0.25))^2.9)
    log10(x+const)
    }
