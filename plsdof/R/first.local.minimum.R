first.local.minimum<-function(x){
    if (length(x)<=2){
        dummy<-which.min(x)
    }
    if (length(x)>2){
        m<-length(x)
        ascending.afterwards<-(x[1:(m-1)]<=x[2:m])
        dummy<-m
        if (sum(ascending.afterwards)>0){
        dummy<-(1:(m-1))[ascending.afterwards]
        dummy<-min(dummy)
        }

    }

return(dummy)
}