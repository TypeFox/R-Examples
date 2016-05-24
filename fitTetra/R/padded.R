padded <-
function(x,maxx) {
  paste(substr("00000",1,nchar(as.character(maxx))-nchar(as.character(x))),x,sep="")
}
