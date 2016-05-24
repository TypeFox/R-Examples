print.npar.CDM <-
structure(function(x, ...)  
{
  output <- x$att.dist
   cat("ACTCD: Asymptotic Classification Theory for Cognitive Diagnosis\n")
  cat("-------------------------------------------\n")
  cat(paste(paste("Analysis starts at", x$starting.time),"\n"))
  cat(paste(paste("Analysis ends at", x$end.time),"\n"))
  cat(paste(paste("based on", x$cluster.method, "cluster algorithm and", x$label.method, "label method"),"\n"))
  cat("-------------------------------------------\n")
  cat("The distribution of estimated attribute patterns\n")
  print(output)
}, export = FALSE, S3class = "npar.CDM", modifiers = "public")
