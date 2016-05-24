`HurstK` <-
function(z){
    y<-z-mean(z)
    S<-cumsum(y)
    R<-(max(S)-min(S))/sqrt(sum(y^2)/length(y))
    log(R)/log(length(y)/2)
}

