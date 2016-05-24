plloss <- function(event, f, Ri) {
    n <- length(event)
    if (length(f) == 1) {
        f <- rep(f, n)
    }
    ef <- exp(f)
    risk <- do.call(c, lapply(X = Ri, FUN = helpfunctionmultistate1, ef = ef))
    lpl <- sum(event * (f - log(risk)))
    return(lpl)}
