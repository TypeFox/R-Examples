print.fdiscd.misclass <-
function(x, dist.print=FALSE, prox.print=FALSE, digits=2, ...)
{
# Misclassification ratio and confusion matrix
alloc <- x$confusion.mat
rsums <- rowSums(alloc)
csums <- colSums(alloc)
wellclass <- diag(alloc)
namesalloc <- names(dimnames(alloc))

alloc <- cbind(alloc, total = rsums, misalloc = round(1 - wellclass/rsums, 3))
alloc <- rbind(alloc, total = c(csums, sum(csums), NA))
names(dimnames(alloc)) <- namesalloc

cat("\n\nmisallocation ratio: ", x$misclassed, "\n\n")

print(alloc, na.print = "", ...)
cat("---------------------------------------------------------------\n")

# Data frame of prior and predicted classes
print(x$classification, digits=3, ...)
if (dist.print)
 {# Matrix of groups-classes distances
  cat("---------------------------------------------------------------\n")
  cat("distances between groups and classes\n"); print(x$distances, digits = digits, ...)}
if (prox.print)
 {# Matrix of groups-classes proximity measures
  cat("---------------------------------------------------------------\n")
  cat("groups-classes proximity index\n"); print(x$proximities, ...)}
return(invisible(x))
}
