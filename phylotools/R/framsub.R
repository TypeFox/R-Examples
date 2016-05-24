#### Function framsub as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

framsub <- function(x, pattern = "-", replacement = "?"){
     result <- rep(NA, 2 * nrow(x))
     dim(result) <- c(nrow(x),2)
     for(i in 1:nrow(x)){
         result[i, 2] <- edgesub(x = x[i, 2], pattern = pattern, replacement = replacement) 
     }
     result[,1] <- x[,1]
     return(result)
}

