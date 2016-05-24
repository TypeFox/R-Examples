transpose <- function(string){
  len <- nchar(string)
  i<-1:(len-1)
  bucket<-c()
  bucket<-append(bucket,unique(as.vector(sapply(string, function(x) paste0(str_sub(string,1,i-1),str_sub(string,i+1,i+1),str_sub(string,i,i),str_sub(string,i+2,nchar(string)))))))
  bucket
}
