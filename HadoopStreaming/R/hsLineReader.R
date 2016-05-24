`hsLineReader` <-
function(file="",chunkSize=-1L, skip=0L, FUN=function(x) cat(x,sep='\n')) {
  if (skip>0) {
    junk = readLines(file, n=skip)
  }
  repeat {
    a = readLines(con=file, n=chunkSize,warn=TRUE)
    if (length(a) ==0) break
    FUN(a)
  }
}

