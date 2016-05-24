`comp.psm` <-
function(cls) {
  
    n <- ncol(cls)
    M <- nrow(cls)

    if(sum(cls %in% (1:n)) < n * M){
        stop("All elements of cls must be integers in 1:nobs")
    }   
    cls.C <- as.vector(t(cls))
    psm <- rep(0,n^2)
    psm <- matrix(.C("comp_psm", as.integer(psm), as.integer(cls.C), as.integer(n), as.integer(M))[[1]]/M,ncol=n) 
    psm
}
