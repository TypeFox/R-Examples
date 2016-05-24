print.summary.simple.surv.sim <-
function(x, ...)
  {
    attr(x,"class") <- "data.frame"
    nn <- c("Number of subjects at risk", "Number of events", "Proportion of subjects with event", 
            "Total time of follow-up", "Time of follow-up (median)", "Density of incidence")
    for (i in 1:length(x))
    {
      cat("\n")
      cat(nn[i], "\n")
      cat("----------------------------\n ")
      print(cbind(x[i]), na.print="", row.names=FALSE,quote=FALSE)
      cat("\n")
    }
  }
