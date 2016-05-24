"rtcmvnorm"<-function(n=1, mean = 0, V = 1, x=0, keep=1, lower = -Inf, upper = Inf){
 
    if(!is.matrix(V)){V<-as.matrix(V)}

    if(length(mean)!=nrow(V) | length(mean)!=ncol(V)){
      stop("V must have same dimensions as mean")
    }
    if(!is.positive.definite(V)){stop("V must be positive definite")}
    if(upper<lower){stop("lower trunaction point has to be less than the upper truncation point")}
    if(keep<0 | keep>nrow(V)){stop(paste("keep must be between 1 and", nrow(V)))}
    if(length(x)!=nrow(V)){stop(paste(" x must have length", nrow(V)))}
    if(lower==-Inf){lower<--1e+35}
    if(upper==Inf){upper<-1e+35}

    rv<-1:n

    output<-.C("rtcmvnormR",
      as.integer(n),
      as.double(mean),
      as.double(x),
      as.double(V),
      as.integer(keep-1),
      as.integer(nrow(V)),	   
      as.double(lower),
      as.double(upper),
      as.double(rv)
    )
    return(output[[9]])
}
