`k_hyd` <-
function(x){
    d<-diff(x)
    res<--c(d,NA)/x
    res[d>=0]<-NA
    res[x==0]<-NA
    return(res)
}

