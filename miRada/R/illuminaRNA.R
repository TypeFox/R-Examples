##  for Illumina microarrays for RNA

print.ilRNA <- function(x,...){
  cat("\nData: ", x$data.name, 
      " (", length(x$gene.names), " RNAs)\n\n", sep = "")
  tmp = rbind(x$mem.x,x$mem.y)
  row.names(tmp) = c("Green","Red")
  cat("\nMeasurement error model parameters:\n ")
  print(tmp)
  cat("\n")
  print(summary(as.data.frame(x[c("x", "y")])), ...)
  invisible(x)
}

