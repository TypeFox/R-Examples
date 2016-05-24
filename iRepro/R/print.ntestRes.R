print.ntestRes <-
function(x, ...){
  cat("Pearson's Chi-squared test for normality of residuals\n")
  t <- as.data.frame(cbind(x$p.value.res1,x$p.value.res2))
  rownames(t) <- "p.value"
  colnames(t) <- c("time = 1", "time = 2")
  print(t)
}

