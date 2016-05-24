"dirnev" <-
function () {

    .K <- getwd()
    .L <- strsplit(.K, "RWork")[[1]]
    .L <- rev(.L)[1]
    NC <- nchar(.L)
    KI <- substring(.L, 2, NC - 2)
    return(KI)
}

