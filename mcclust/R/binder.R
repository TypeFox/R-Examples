`binder` <-
function(cls,psm){
    if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
     stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
    
    if(is.vector(cls)) cls <- t(cls)
    apply(cls,1, function(x) sum(abs(cltoSim(x)-psm))/2)
}

