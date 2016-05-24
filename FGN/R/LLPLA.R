`LLPLA` <-
function(alpha, z){
    if (alpha<0.01 || alpha> 1.99)
        return(-10^10)
    r<-acvfPLA(alpha, length(z)-1)
    DLLoglikelihood(r, z)
}
