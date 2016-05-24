is.vlabeled <-
function(x){
  
  l <- memisc::labels(x)
  
  if (is.null(l)) {
    vl <- 999999
  } else {
  vl <- l@values
  }
  x %in% vl
}
