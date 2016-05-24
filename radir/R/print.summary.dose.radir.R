print.summary.dose.radir <-
function(x, ...)
{
  attr(x,"class") <- "data.frame"
  nn <- c("Mode", "Expected value", "Standard Dev.", "95% CI")
  for (i in 1:4)
  {
    cat("\n")
    cat(nn[i], "\n")
    cat("----------------------------\n ")
    write.table(x[i], quote=FALSE, na="NA", row.names=FALSE, col.names=FALSE)
    cat("\n")
  }
}
