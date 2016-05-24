`pear` <-
function(cls,psm){
    if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
     stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
    
    if(is.vector(cls)) cls <- t(cls)
    
    compIpi <- function(cl){
       mat <- cltoSim(cl)*psm
       sum(mat[lower.tri(mat)])
    }
    
    n <- dim(psm)[1]
    no2 <- choose(n,2)
    sumpij <- sum(psm[lower.tri(psm)])
    sumIij <- apply(cls,1, function(x) sum(choose(table(x),2)))
    sumIpi <- apply(cls,1,compIpi)
    correc <- (sumIij*sumpij)/no2
    (sumIpi - correc)/ (0.5*(sumpij+sumIij)-correc)
}

