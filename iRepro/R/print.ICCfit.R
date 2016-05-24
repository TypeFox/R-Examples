print.ICCfit <-
function(x, ...){
  t <- as.data.frame(x$icc)
  rownames(t) <- "ICC"
  colnames(t) <- ""
  print(t)
}

