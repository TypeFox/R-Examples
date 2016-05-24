#### Function formatXY as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


formatXY <- function(x){
    ncharX <- substring(x, 2, regexpr("Y", x)-1)
    ncharY <- substring(x, nchar(ncharX)+3, nchar(x))
    resX <- c()
    for(i in 1:length(ncharX)){
        n0x <- paste( rep(rep(0, length(ncharX[i])), 
                     times = (max(nchar(ncharX)) - nchar(ncharX[i]) + 1)),
                     collapse = "", sep = "")
        resX[i] <- paste("X", substring(n0x, 2, nchar(n0x)), ncharX[i], collapse = "", sep = "")
    }
    resY <- c()
    for(i in 1:length(ncharY)){
        n0x <- paste( rep(rep(0, length(ncharY[i])), 
                     times = (max(nchar(ncharY)) - nchar(ncharY[i]) + 1)),
                     collapse = "", sep = "")
        resY[i] <- paste("Y", substring(n0x, 2, nchar(n0x)), ncharY[i], collapse = "", sep = "")
    }
    res <- paste(resX, resY, sep = "")
    return(res)
}

