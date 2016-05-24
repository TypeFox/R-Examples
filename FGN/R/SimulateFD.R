`SimulateFD` <-
function(n, d){
    stopifnot(d < 0.5, n>1)
    r<-acvfFD(d, n-1)
    z<-try(DHSimulate(n, r), silent=TRUE)
    if(class(z)=="try-error") z<-DLSimulate(n, r)
    z
}
