is.numeric_data.frame=function(x){
  if (is.data.frame(x) && all(sapply(x,base::is.numeric)))
    return (T)
  return (F)
}

is.numeric=function(x){
  if (base::is.numeric(x))
    return (T)
  if (is.data.frame(x) && all(sapply(x,base::is.numeric)))
    return (T)
  return (F)
}

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}