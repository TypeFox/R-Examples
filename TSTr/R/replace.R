replace <- function(string, useUpper){
  len <- nchar(string)
  i<-1:len
  bucket<-c()
  if(useUpper == T){
    letter<-c(letters,LETTERS,c(0,1,2,3,4,5,6,7,8,9))
  }
  else{
    letter<-c(letters,c(0,1,2,3,4,5,6,7,8,9))
  }
  for(letra in letter){
    bucket<-append(bucket,as.vector(sapply(string, function(x) paste0(str_sub(string,1,i-1),letra,str_sub(string,i+1,nchar(string))))))
  }
  bucket<-unique(bucket)
  bucket
}
