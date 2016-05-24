obsprof.wprof <-
function(prof) {
    res <- list()
    wh <- which(prof$freq > 0)
    res <- list(profiles=prof$profiles[wh,], freq=prof$freq[wh])
    class(res) <- "wprof"
    return(res)
}
