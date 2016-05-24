print.mcp2 <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  df <- x$comp
  gnames1 <- x$fnames[x$comp[,1]]
  gnames2 <- x$fnames[x$comp[,2]]
  rownames(df) <- paste(gnames1, gnames2, sep = " vs. ")
  colnames(df)[1:2] <- paste("Group", 1:2)
  
  if (colnames(df)[7] == "p.crit") {
      sig <- ifelse(x$comp[,6] < x$comp[,7], TRUE, FALSE)
      df <- round(df, 5)
      df <- data.frame(df, sig)
      print(df[,-c(1:2)])
    } else {
      print(as.data.frame(round(df[,-c(1:2)], 5)))
  }
  
  cat("\n")
}
