### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2011-02-19 21:43:29 umusus>
### jyear

jyear <- function(x,shift=FALSE,withAD=FALSE,ascii=FALSE){
  if(!is.numeric(x)) x <- as.numeric(x)
  if(x < 1868) stop(paste(x," is not year or too ancient"))
  if(!(x == 1912 || x == 1926 || x == 1989) && shift) {
    stop(paste("'shift' is not for year:",x))
  }
  gengoKanji <- jdate$gengo
  gengoAscii <- c("M", "T", "S", "H")
  if (ascii) {
    gengo <- gengoAscii
    suffix <- NULL
  }
  else {
    gengo <- gengoKanji
    suffix <- jdate$ymd[1]
  }
  i <- ifelse(shift,-1,0)
  if (x < 1912) {
    y <- x - 1867
    g <- gengo[1]
  }
  else if (x < 1926) {
    y <- x - 1911 + ifelse(shift,44,0)
    g <- gengo[2 + i]
  }
  else if (x < 1989) {
    y <- x - 1925 + ifelse(shift,14,0)
    g <- gengo[3 + i]
  }
  else {
    y <- x - 1988 + ifelse(shift,63,0)
    g <- gengo[4 + i]
  }
  if(withAD) res <- paste(g,y,"(",x,")",suffix,sep='')
  else res <- paste(g,y,suffix,sep='')
  return(res)
}

