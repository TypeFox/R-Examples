subNonStandardCharacters <- function(x,
   standardCharacters=c(letters, LETTERS, ' ','.', '?', '!', 
      ',', 0:9, '/', '*', '$', '%', '\"', "\'", '-', '+', '&', '_', ';', 
      '(', ')', '[', ']', '\n'),
   replacement='_',
   gsubList=list(list(pattern='\\\\\\\\|\\\\', replacement='\"')),
   ... ) {
##
## 1.  length(x)<1?
##
  nx <- length(x)
  if(nx<1)return(integer(0))
##
## 2.  gsubList
##
  xo <- x
  ng <- length(gsubList)
  for(ig in seq(length=ng)){
      gsLi <- gsubList[[ig]]
      if(!is.list(gsLi)){
          print(deparse(substitute(gsubList)))
          cat(ig, '\n')
          print(gsLi)
          stop('gsubList[[', ig, ']] is NOT a list')
      }
      xo <- gsub(gsLi$pattern, gsLi$replacement, xo)
  }
##
## 3.  stringi (formerly iconv)
##
#  xo <- iconv(xo, "", "ASCII//TRANSLIT")
#  returns NAs on Mac OS X 10.10.1 2015.01.25 
  Encoding(xo)[Encoding(xo)=='unknown'] <- 'latin1'
  xi <- stringi::stri_trans_general(xo, "Latin-ASCII")
##
## 4.  x. <- strsplit(x, "", ...)
##
  x. <- strsplit(xi, "", ...)
##
## 5.  check each and modify as needed
##
  for(ix in seq(length=nx)){
      gi <- which(!(x.[[ix]] %in% standardCharacters))
      if(length(gi)>0){
        if(!is.na(xi[ix])){
          gi. <- range(gi)
          ni <- length(x.[[ix]])
          xi1 <- paste(x.[[ix]][seq(length=gi.[1]-1)], collapse='')
          xie <- paste(x.[[ix]][seq(from=gi.[2]+1, length=ni-gi.[2])],
                       collapse='')
          xi[ix] <- paste(xi1, replacement, xie, sep="")
        }
      }
  }
##
## 5.  Done
##
  xi
}
