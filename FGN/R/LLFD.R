`LLFD` <-
function(d, z){
    if (abs(d)>0.49) 
        return(-10^10)
    r<-acvfFD(d, length(z)-1)
    DLLoglikelihood(r, z)
}

