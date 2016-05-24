"rtnorm"<-function(n=1, mean = 0, sd = 1, lower = -Inf, upper = Inf){
 
    if(any(lower==-Inf)){lower[which(lower==-Inf)]<-c(-1e+35)}
    if(any(upper==Inf)){upper[which(upper==Inf)]<-1e+35}
    if(any(sd<=0)){stop("standard deviation must be positive")}
    if(length(mean)==1){
      mean<-rep(mean, n)
    }
    if(length(sd)==1){
      sd<-rep(sd, n)
    }
    if(length(lower)==1){
      lower<-rep(lower, n)
    }
    if(length(upper)==1){
      upper<-rep(upper, n)
    }
    if(any(c(length(mean), length(sd), length(lower), length(upper))!=n)){stop("inputs must be the same length or of length 1")}
    if(any(upper<lower)){stop("lower trunaction point has to be less than the upper truncation point")}
    rv<-1:n

    output<-.C("rtnormR",
      as.integer(n),
      as.double(mean),
      as.double(sd),			   
      as.double(lower),
      as.double(upper),
      as.double(rv)
    )
    return(output[[6]])
}
