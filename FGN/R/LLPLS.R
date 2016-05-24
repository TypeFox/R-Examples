`LLPLS` <-
function(alpha, z){
    if (alpha<0.001 || alpha> 0.999)
        return(-10^10)
    r<-acvfPLS(alpha, length(z)-1)
    DLLoglikelihood(r, z)
}
