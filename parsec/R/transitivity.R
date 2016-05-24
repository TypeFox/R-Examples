transitivity <-
function(m) {
    r <- nrow(m)
    c <- ncol(m)
    if(r!=c) return(FALSE)
    all(((m %*% m) > 0) <= m)
}
