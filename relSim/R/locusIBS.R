locusIBS = function(profMat){

    ## profIbs = function(prof1, prof2){
    ##     res = 0

    ##     a1 = prof1[1]
    ##     a2 = prof1[2]
    ##     b1 = prof2[1]
    ##     b2 = prof2[2]

    ##     if(a1 == b1 && a2 == b2){
    ##         res = 2
    ##     }else if((a1 == b1) || (a2 == b1) || (a1 == b2) || (a2 == b2)){
    ##         res = 1
    ##     }else{
    ##         res = 0
    ##     }

    ##     return(res)
    ## }

    N = nrow(profMat)
    nc = ncol(profMat)

    if(nc != 4)
        stop("Wrong dimensions")

    p = as.vector(t(profMat))
    return(.locusIBS(p, N))
}
