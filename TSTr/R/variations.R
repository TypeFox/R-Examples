variations <- function(vector, useUpper){
  bucket<-c()
  for(word in vector){
    x <- append(deletion1(word),transpose(word))
    x <- append(x, replace(word, useUpper))
    x <- append(x, insertion(word, useUpper))
    x <- unique(x)
    bucket<-append(bucket,x)
  }
  bucket
}
