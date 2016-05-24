'make.v'<-function(n,r,sig2){
     out<-matrix(r,n,n)
     for (i in 1:n){
          out[i,i]<-1
     }
    out<-sig2*out
    return(out)
}
