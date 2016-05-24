print.wilcoxh <- function(x, ...)
{
  #print method for objects of class wilcoxh (from WilcoxH)
  res.table <- as.data.frame(rbind(x$SH.res, x$WH.res))
  row.names(res.table) <- c("Steiger-Hakstian","Wilcox H")
  colnames(res.table) <- c("X^2","df","p-value")

  cat("\nResults for 0-correlation tests: \n")
  print(res.table)
  cat("\n")
}
