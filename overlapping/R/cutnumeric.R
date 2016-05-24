cutnumeric <-
function(x,n=1000) {
  ## x = a numeric vector
  xclass <- cut(x,seq(min(x),max(x),length=n),include.lowest=TRUE)
  xnumeric <- NULL
  for (l in xclass) {
    h <- strsplit(l,",")[[1]]
    h <- gsub("\\[","",gsub("\\]","",gsub("\\(","",h)))
    xnumeric <- c(xnumeric,mean(as.numeric(h)))
  } 
  return(xnumeric)
}
