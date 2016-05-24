`LLFGN` <-
function(H, z){
    if (H<0.01 || H>0.99) 
        return(-10^10)
    DLLoglikelihood(acvfFGN(H, length(z)-1), z)
}

