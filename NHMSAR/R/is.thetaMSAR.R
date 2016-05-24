is.thetaMSAR <-
function(x) {
    if (!is.list(x)) {
        return(FALSE)
    }
    if (length(x) < 4) {
        FALSE
    } else {
        TRUE
    }
}
