"dcmvnorm"<-function(x, mean = 0, V = 1, keep=1, cond=(1:length(x))[-keep], log=FALSE){
 
    if(!is.matrix(V)){V<-as.matrix(V)}

    n<-nrow(V)
    if(ncol(V)!=n){stop("V must be square")}

    if(length(mean)!=nrow(V) | length(mean)!=ncol(V)){
      stop("V must have same dimensions as mean")
    }
    if(length(x)!=nrow(V) | length(mean)!=ncol(V)){
      stop("V must have same dimensions as mean")
    }
    if(!is.positive.definite(V)){stop("V must be positive definite")}
    if(any(keep<0 | keep>nrow(V))){stop(paste("keep must be between 1 and", nrow(V)))}
    if(any(cond<0 | cond>nrow(V))){stop(paste("cond must be between 1 and", nrow(V)))}
    if(any(keep%in%cond)){stop("keep and cond should be distinct")}   

    output<-.C("dcmvnormR",
      as.integer(n),
      as.double(x),
      as.double(mean),
      as.double(V),
      as.integer(keep-1),
      as.integer(cond-1),	   
      as.integer(length(keep)),
      as.integer(length(cond)),
      as.double(1)
    )
    if(!log){
      output[[9]]<-exp(output[[9]])
    }
    return(output[[9]])
}
