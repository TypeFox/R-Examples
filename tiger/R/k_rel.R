`k_rel` <-
function(x,y){
    rel <- k_hyd(x) / k_hyd(y)
    toRet<-mean(rel,na.rm=TRUE)
    return(toRet)
}

