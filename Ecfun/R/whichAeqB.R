whichAeqB <- function(A, B, errNoMatch='no match',
                      err2Match='more than one match'){
    ab <- which(A %in% B)
    nab <- length(ab)
    if(nab<1){
        stop(errNoMatch)
    }
    if(nab>1){
        stop(err2Match)
    }
    ab
}
