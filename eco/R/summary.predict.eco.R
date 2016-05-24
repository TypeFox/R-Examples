summary.predict.eco <- function(object, CI=c(2.5, 97.5), ...) {

  if (any(CI < 0) || any(CI > 100))
    stop("Invalid input for CI")
  n.draws <- nrow(object)
  n.var <- ncol(object)
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
			paste(max(CI), "%", sep=" "))
  W.table <- matrix(NA, ncol=length(table.names), nrow=n.var)
  for (i in 1:n.var)
    W.table[i,] <- cbind(mean(object[,i]), sd(object[,i]), 
                         quantile(object[,i], min(CI)/100), 
                         quantile(object[,i], max(CI)/100))
  colnames(W.table) <- table.names
  rownames(W.table) <- colnames(object)

  res <- list(W.table = W.table, n.draws = n.draws)
  class(res) <- "summary.predict.eco"
  return(res)
}
