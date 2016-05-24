print.summary.crisk.data.sim <-
function(x, ...)
{
    attr(x,"class") <- "data.frame"
    nn <- c("", "Number of subjects at risk", "Number of events", "Density of incidence")
    for (i in 2:length(x))
    {
      cat("\n")
      cat(nn[i], "\n")
      cat("----------------------------\n ")
      print(cbind(x[1],x[i]), na.print="", row.names=FALSE,quote=FALSE)
      cat("\n")
    }
}
