grepNonStandardCharacters <- function(x, value=FALSE,
   standardCharacters=c(letters, LETTERS, ' ','.', ',', 0:9,
       '\"', "\'", '-', '_', '(', ')', '[', ']', '\n'),
   ...) {
##
## 1.  split single characters
##
  if(is.factor(x)){
      x <- as.character(x)
  }
  if(!is.character(x)){
      stop('x is not character;  is ', class(x))
  }
  x. <- strsplit(x, '', ...)
##
## 2.  find !standardCharacters
##
  x.in.sC <- function(y){
      x0 <- which(!(y %in% standardCharacters))
      if(length(x0)<1) 0 else x0[1]
  }
  x2 <- sapply(x., x.in.sC)
##
## 3.  Return indices or value
##
  x1 <- which(x2>0)
  if(value){
      return(x[x1])
  } else return(x1)
}
