print.cd.cluster <-
structure(function(x, ...)  
{
  output1 <- cbind(c(1:length(x$size)),x$size)
  colnames(output1) <- c("clusters #","freq")
  cat("-------------------------------------------\n")
  cat("Cluster analysis for ACTCD\n")
  cat(paste(paste("based on", x$cluster.method, "algorithm"),"\n"))
  cat("-------------------------------------------\n")
  cat("Number of examinees within each cluster:\n")
  print(output1)

}, export = FALSE, S3class = "cd.cluster", modifiers = "public")
