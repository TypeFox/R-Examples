
.onAttach <- function(lib, pkg) {
    packageStartupMessage(paste(" Package", sQuote("exactRankTests"), 
        "is no longer under development.\n",
        "Please consider using package", sQuote("coin"), "instead.\n"))
}

cpermdist2 <- function(N, totsum, sc, scores, log) {
    .Call("cpermdist2", as.integer(N),
          as.integer(totsum), as.integer(sc),
          as.integer(scores), as.logical(FALSE),
          PACKAGE="exactRankTests")
}
