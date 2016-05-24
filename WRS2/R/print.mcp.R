print.mcp <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  if (length(x$effects[[1]][[1]]) > 1) facA <- cbind(x$effects[[1]][[1]], x$effects[[1]][[2]], x$effects[[1]][[3]]) else facA <- c(x$effects[[1]][[1]], x$effects[[1]][[2]], x$effects[[1]][[3]])
  if (length(x$effects[[2]][[1]]) > 1) facB <- cbind(x$effects[[2]][[1]], x$effects[[2]][[2]], x$effects[[2]][[3]]) else facB <- c(x$effects[[2]][[1]], x$effects[[2]][[2]], x$effects[[2]][[3]])
  if (length(x$effects[[3]][[1]]) > 1) facAB <- cbind(x$effects[[3]][[1]], x$effects[[3]][[2]], x$effects[[3]][[3]]) else facAB <- c(x$effects[[3]][[1]], x$effects[[3]][[2]], x$effects[[3]][[3]])
  df <- rbind(facA, facB, facAB)
  rownames(df) <- colnames(x$contrasts)
  colnames(df)[4] <- "p-value"
  print(as.data.frame(round(df, 5)))
  
  if(exists("x$alpha.crit")) cat("\nThe critical alpha level is ", x$alpha.crit, ".", sep = "")
  cat("\n")
}
