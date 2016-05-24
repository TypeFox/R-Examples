#### Function dimension as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

dimension <-
function(x, unique = FALSE, sort = FALSE){
    ncharX <- substring(x, 2, regexpr("Y", x)-1)
    ncharY <- substring(x, nchar(ncharX)+3, nchar(x))
    if(unique){
        ncharX = unique(ncharX)
        ncharY = unique(ncharY)
    }
    if(sort){
        ncharX = sort(as.numeric(ncharX))
        ncharY = sort(as.numeric(ncharY))
    }
    res <- list(ncharX, ncharY)
    return(res)
}

