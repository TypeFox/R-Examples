normalize <-
function(v,w=NULL){
    norm.v<-sqrt(sum(v^2))
    v<-v/norm.v
    if (is.null(w)==TRUE){
        return(v)
    }
    if (is.null(w)==FALSE){
        w<-w/norm.v
        return(list(v=v,w=w))
    }
}

