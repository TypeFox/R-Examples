print.mcp1 <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  df <- x$comp
  gnames1 <- x$fnames[x$comp[,1]]
  gnames2 <- x$fnames[x$comp[,2]]
  rownames(df) <- paste(gnames1, gnames2, sep = " vs. ")
  
  if (ncol(df) > 6) {
    if (colnames(df)[7] == "crit") {
      sig <- ifelse(x$comp[,6] > x$comp[,7], TRUE, FALSE)
      df <- round(df, 5)
      df <- data.frame(df, sig)
      colnames(df)[6] <- "test"
      print(df[,-c(1:2)])
    } else {
      print(as.data.frame(round(df[,-c(1:2)], 5)))
    }
  }
  else {
      print(as.data.frame(round(df[,-c(1:2)], 5)))
    }
  
  cat("\n")
}
