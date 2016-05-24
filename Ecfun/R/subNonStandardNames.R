subNonStandardNames <- function(x,
   standardCharacters=c(letters, LETTERS, ' ','.', '?', '!', 
      ',', 0:9, '/', '*', '$', '%', '\"', "\'", '-', '+', '&', '_', ';', 
      '(', ')', '[', ']', '\n'),
   replacement='_',
   gsubList=list(list(pattern='\\\\\\\\|\\\\', replacement='\"')),
   removeSecondLine=TRUE, nonStandardNames=Ecdat::nonEnglishNames, 
   namesNotFound="attr.replacement", ...) {
##
## 1.  removeSecondLine
##
  x0 <- x
  if(is.data.frame(x0))x <- as.matrix(x0)
  if(length(x)<1){
    return(NULL)
  }
  if(removeSecondLine){
    nch0 <- (nchar(x)<1)
    X2. <- strsplit(x, '\n')
    x2 <- sapply(X2., function(xx){
        xx2 <- xx[-1]
        paste(xx2, collapse='\n')
    } ) 
    no2 <- sapply(X2., length)
    x2[no2<2] <- NA
    x <- sapply(X2., '[', 1)
    x[nch0] <- '' 
  } else x2 <- NULL
##
## 2.  x. <- subNonStandardCharacters(x, ...)
##
  x. <- subNonStandardCharacters(x, standardCharacters, replacement,
                                 gsubList, ...)
##
## 3.  loop over rows of nonStandardNames
##
  nSN <- nrow(nonStandardNames)
  for(iSN in seq(length=nSN)){
      x. <- gsub(nonStandardNames[iSN, 1], nonStandardNames[iSN, 2],
                 x.)
  }
##
## 4.  Eliminate leading and trailing blanks
##
#  if(require(tis)){
      x. <- tis::stripBlanks(x.)
#  } else {
#      warning('need stripBlanks{tis} to delete leading and trailing',
#              ' blanks;  not available')
#  }
##
## 5.  namesNotFound?
##
  if('attr.replacement' %in% namesNotFound){
    nNF <- grep(replacement, x.) 
  } else if('attr.notFound' %in% namesNotFound){
    nNF <- which(x0 != x.) 
  } else nNF <- integer(0)
##
## 6.  Reformat as matrix or data.frame 
##
  x1 <- x. 
  if(is.matrix(x0)){
    attributes(x.) <- attributes(x0)
  } else if(is.data.frame(x0)){
    dim(x.) <- dim(x0)
    colnames(x.) <- colnames(x0)
    x. <- as.data.frame(x., stringsAsFactors=FALSE)
  }
##
## 7.  attr(x., 'namesNotFound') 
##
  if(length(nNF)>0){
    attr(x., 'namesNotFound') <- x1[nNF]
    if('print' %in% namesNotFound){
      cat('Non-standard names not in nonStandardNames: ', 
          x1[nNF], '\n')
    }
  }
##
## 8.  Done
##
  nchx2 <- nchar(x2)
  nchx2[is.na(x2)] <- 0 
  if(any(nchx2>0))attr(x., 'secondLine') <- x2
  x.
}
