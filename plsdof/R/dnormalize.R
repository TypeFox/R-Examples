dnormalize <-
function(v,dv){
    n<-length(v)
    vn<-normalize(v)
    dn=(1/sqrt(sum(v^2)))*(diag(n)- vn%*%t(vn))%*%dv
    return(dn)
}

