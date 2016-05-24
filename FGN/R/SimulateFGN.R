`SimulateFGN` <-
function(n, H){
    stopifnot( H>0, H<1, n>1)
    r<-acvfFGN(H, n-1)
    if (n>=50 && H<0.84)
        z<-DHSimulate(n, r)
    else
        z<-DLSimulate(n, r)
    z
}
