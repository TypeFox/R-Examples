deletion2 <- function(string){
  len <- nchar(string)
  i<-1:len
  j<-1:(len-1)
  bucket<-c()
  bucket<-append(bucket,string)
  bucket<-append(bucket,unique(as.vector(sapply(string, function(x) paste0(str_sub(string,1,i-1),str_sub(string,i+1,nchar(string)))))))
  bucket<-append(bucket,unique(as.vector(sapply(as.vector(sapply(string, function(x) paste0(str_sub(x,1,i-1),str_sub(x,i+1,nchar(x))))), function(x) paste0(str_sub(x,1,j-1),str_sub(x,j+1,nchar(x)))))))
  bucket
}
