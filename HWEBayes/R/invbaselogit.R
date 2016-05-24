invbaselogit <-
function(baselogit) {
    p <- c(exp(baselogit), 1)
    list(probs = p/sum(p))
}

