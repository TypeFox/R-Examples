"julian2date" <-
  function(x){
    orig <- as.Date(attributes(x)[[1]])
    jorig <- as.numeric(orig)
    seqdates <- seq(from=orig,to=orig+max(x, na.rm=TRUE),by=1)
    seqjulian <- seq(from=jorig,to=jorig+max(x, na.rm=TRUE),by=1)
    seqdates[x+1]
}
