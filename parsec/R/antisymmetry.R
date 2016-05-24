antisymmetry <-
function(m) {
    r <- nrow(m)
    c <- ncol(m)
    if(r!=c) return(FALSE)
    all((m+t(m) - diag(diag(m))) <= 1)
}
