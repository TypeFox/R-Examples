dztp <-
function(y,lambda){
    n<-length(y)
    out<-rep(0,n)
    out[y>=1]<-dpois(y[y>=1],lambda)/(1-dpois(0,lambda))    
return(out)
}
