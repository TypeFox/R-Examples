print.labeling <-
structure(function(x, ...)  
{
  output2 <- x$att.dist
    
  cat("-------------------------------------------\n")
  cat("labeling for ACTCD\n")
  cat(paste(paste("based on", x$label.method, "label method"),"\n"))
  cat("-------------------------------------------\n")
  cat("The distribution of attribute patterns:\n")
  print(output2)
 
}, export = FALSE, S3class = "labeling", modifiers = "public")
