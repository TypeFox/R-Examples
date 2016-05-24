doQmap <- function(x,fobj,...){
  cc <- class(fobj)
  ffun <- substring(cc,4,nchar(cc))
  ffun <- paste("do",ffun,sep="")
  test <- sapply(ffun,exists,mode="function") 
  if(all(test)){
    ffun <- match.fun(ffun)
  } else {
    stop("doQmap not defined for class(fobj) ==",
         class(fobj))
  }
  ffun(x,fobj,...)
}
