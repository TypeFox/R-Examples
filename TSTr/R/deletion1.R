deletion1 <- function(string){
  len <- nchar(string)
  i<-1:len
  bucket<-c()
  bucket<-append(bucket,string)
  bucket<-append(bucket,unique(as.vector(sapply(string, function(x) paste0(str_sub(string,1,i-1),str_sub(string,i+1,nchar(string)))))))
  bucket
}
