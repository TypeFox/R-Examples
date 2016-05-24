adj.input <-
function(x,threshold=.5) {
  
  if(memisc::measurement(x) %in% c("nominal","ordinal")){
    
    nlabeled <- sum(is.vlabeled(x),na.rm=TRUE)
    if(nlabeled < threshold*length(x))
      memisc::measurement(x) <- "interval"
  }
  x
}
