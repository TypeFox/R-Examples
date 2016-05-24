g.test <- function(y, correct = FALSE, pi.null = NULL){
  c <- chisq.test(y, correct = correct)
  E <- ifelse(is.null(pi.null), c$expected, pi.null * sum(y)) 
  temp <- log(y/E)
  temp <- ifelse(temp==-Inf,0, temp)

  G <- sum(y * temp) * 2
  Df <- c$parameter
  p.val <- pchisq(G, df = Df, lower.tail=FALSE)
  res <- list()
  res$head <- "Contingency table likelihood ratio test"
  res$summary <- data.frame(G = G, P = p.val)
  names(res$summary) <- c("G.statistic", "P-value")
  row.names(res$summary) <- ""
  class(res) <- "gtest"

  #---------------------------------#

cat("\n")
cat(res$head,"\n")
cat("\n")
rq<-structure(res$summary)
print(rq, digits = max(3, getOption("digits")))
cat("\n")
invisible(res)
}

