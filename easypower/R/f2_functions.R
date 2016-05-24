# Funtions to convert eta.sq to f2
# Converts strings to numeric then calculates the f2 value
#
# Multiple groups function
f2.vector <- function(eta.sq) {
    es.vector <- eta.sq
    for(i in 1:length(es.vector)) {
        if(es.vector[i] == "small") {
            es.vector[i] <- 0.01
        }
        else if(es.vector[i] == "med") {
            es.vector[i] <- 0.06
        }
        else if(es.vector[i] == "large") {
            es.vector[i] <- 0.14
        }
    }
    es.vector <- as.numeric(es.vector)
    out.vector <- sapply(1:length(es.vector), FUN = function(x) {
        es.vector[x]/(1 - es.vector[x])
    }
    )
    return(out.vector)
}

# One-way function
f.oneway <- function(effect.size) {
    if(effect.size == "small") {
        effect.size <- 0.01
    }
    else if(effect.size == "med") {
        effect.size <- 0.06
    }
    else if(effect.size == "large") {
        effect.size <- 0.14
    }

    convert.to.f <- sqrt(effect.size/(1 - effect.size))
    return(convert.to.f)
}
