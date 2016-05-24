file.copy2 <- function(from, to) {
##
## 1.  length(from)
##
  nf <- length(from) 
  if(nf<1)stop("no file to copy")
  if(nf>1)
    stop("'file.copy2' can copy only one file per call;  ",
         "length(from) = ", nf) 
##
## 2.  file.exists(from)?  
##
  if(!file.exists(from))return(FALSE) 
##
## 3.  Check 'to' 
##
  {
    if(missing(to)) to <- from
    else {
      nt <- length(to)
      if(nt<1) to <- from
      else
        if(nt>1)
          stop("'file.copy2' can make only one copy;  ",
               "length(to) = ", nt)
    }
  }
##
## 4.  file.exists(to)?  
##
  if(file.exists(to)){
    Dir <- dir(dirname(to))
    To <- basename(to) 
    nch <- nchar(To)
    DIR <- toupper(Dir) 
    to2 <- (toupper(To) == substring(DIR, 1, nch))
    nto2 <- sum(to2)
    To2 <- paste(To, 1:(nto2+1), sep="")
    Exists <- (toupper(To2) %in% DIR)
    to <- To2[!Exists][1] 
  }
##
## 5.  file.copy   
##
  file.copy(from, to) 
}
