print.fdiscd.predict <-
function(x, dist.print=TRUE, prox.print=FALSE, digits=2, ...)
{
if (!is.null(x$misclassed))
 {# Misclassification ratio and confusion matrix (if they are available in x)
  alloc <- x$confusion.mat
  rsums <- rowSums(alloc)
  csums <- colSums(alloc)
  wellclass <- diag(alloc)
  namesalloc <- names(dimnames(alloc))
  
  alloc <- cbind(alloc, total = rsums, misalloc = round(1 - wellclass/rsums, 3))
  alloc <- rbind(alloc, total = c(csums, sum(csums), NA))
  names(dimnames(alloc)) <- namesalloc
  
  cat("misallocation ratio (computed on the data for which the class is known))\n")
  cat(x$misclassed, "\n\n")
  cat("\n")
  
  print(alloc, na.print = "", ...)
  cat("---------------------------------------------------------------\n")}
# Data frame of prior and predicted classes
print(x$prediction)
if (dist.print)
 {# Matrix of groups-classes distances
  cat("---------------------------------------------------------------\n")
  cat("distances between groups and classes\n"); 
  print(x$distances, digits = digits, ...)}
if (prox.print)
 {# Matrix of groups-classes proximity measures
  cat("---------------------------------------------------------------\n")
  cat("groups-classes proximity index\n"); 
  print(x$proximities, ...)}
return(invisible(x))
}
