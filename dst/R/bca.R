bca<-function(m,f,n){
# State a basic chance assignment over subsets of a frame of discernement
  if((abs(sum(m)-1)>0.000001) | (length(f[,1])!=length(m)) | (length(f[1,])!=length(n))){ 
    stop("error in input arguments: check your input data ") 
    }
  else{
    z<-cbind(m,f)
    colnames(z)<-c("mass",n)
    y<-list(DempsterRule=z,con=0)
    return(y)
  }
 }