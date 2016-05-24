## dt2strO <- function(dt,verbose=FALSE) {
##   hr <- dt %/% 3600; rest <- dt %% 3600
##   mi <- rest %/% 60; rest <- rest %%60
##   se <- rest
##   if (verbose) paste(hr,"hours",mi,"minutes",se,"seconds",sep=" ")
##   else  paste(hr, mi, se, sep = ":")
## }

dt2str <- function(dt,dec=0,verbose=FALSE) {
  hr <- dt %/% 3600; rest <- dt %% 3600
  mi <- rest %/% 60; rest <- rest %%60
  se <- rest
  if (verbose) paste(hr,"hours",mi,"minutes", sprintf(paste("%02.",dec,"f",sep=""), se),"seconds",sep=" ")
  else sprintf(paste("%02d:%02d:%02.",dec,"f",sep=""), hr, mi, se)
}
