entropy.Renyi <-function(wind,q){
    
    H<-apply(wind,2,function(y){
    
        out<-.C("entropyRenyi",
            prob=as.double(y),
            lengthprob=as.double(length(y)),
            q=as.double(q),
            H=double(1))
         out$H
    })
    return(H)
}
