
#########################################################################
#      Realized multipower variation
#########################################################################

setGeneric("mpv",function(yuima,r=2,normalize=TRUE)standardGeneric("mpv"))

setMethod("mpv",signature(yuima="yuima"),
          function(yuima,r=2,normalize=TRUE)mpv(yuima@data,r,normalize))

setMethod("mpv",signature(yuima="yuima.data"),
          function(yuima,r=2,normalize=TRUE){
  
  data <- get.zoo.data(yuima)
  d.size <- length(data)
  
  result <- double(d.size)
  
  if(is.numeric(r)){
    
    for(d in 1:d.size){
      
      X <- as.numeric(data[[d]])
      idt <- which(is.na(X))
      if(length(idt>0)){
        X <- X[-idt]
      }
      if(length(X)<2) {
        stop("length of data (w/o NA) must be more than 1")
      }
      
      abs.diffX <- abs(diff(X))
      tmp <- rollapplyr(abs.diffX,length(r),FUN=function(x)prod(x^r))
      result[d] <- length(abs.diffX)^(sum(r)/2-1)*sum(tmp)
      
    }
    
    if(normalize){
      result <- result/prod(2^(r/2)*gamma((r+1)/2)/gamma(1/2))
    }
    
  }else{
    for(d in 1:d.size){
      
      X <- as.numeric(data[[d]])
      idt <- which(is.na(X))
      if(length(idt>0)){
        X <- X[-idt]
      }
      if(length(X)<2) {
        stop("length of data (w/o NA) must be more than 1")
      }
      
      abs.diffX <- abs(diff(X))
      tmp <- rollapplyr(abs.diffX,length(r[[d]]),FUN=function(x)prod(x^r[[d]]))
      result[d] <- length(abs.diffX)^(sum(r[[d]])/2-1)*sum(tmp)
      
      if(normalize){
        result[d] <- result[d]/prod(2^(r[[d]]/2)*gamma((r[[d]]+1)/2)/gamma(1/2))
      }
    }
  }
  
  return(result)
})