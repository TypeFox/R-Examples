"szelet" <-
function (BE=.QKEP,MIT=c("Z"==substring(names(.QND),1,1)),HOL=c(1,1)) {
  assign("..SZELET",BE,pos=1)
  if(sum(MIT)!=length(HOL)) {
    return(NULL)
  }

  M <- length(MIT)
  KI <- "..SZELET["
  hol <- 1
  for (m in 1:M) {
    if(MIT[m]) {
      mit <- as.character(HOL[hol])
      hol <- hol + 1
    }
    else {
      mit <- ""
    }
    KI <- paste(sep="",KI,mit,",")
  }
  KI <- substring(KI,1,nchar(KI)-1)
  KI <- paste(sep="",KI,"]")
  eva(KI)
}

