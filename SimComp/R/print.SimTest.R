print.SimTest <-
function(x,digits=4,...) {


cat("", "\n")
cat("Test for", x$test.class, "of means of multiple endpoints", "\n")
cat("Assumption: ")
  if (x$covar.equal==TRUE) cat("Homogeneous ") else cat("Heterogeneous ")
  cat("covariance matrices for the groups", "\n")
cat("Alternative hypotheses: True", x$test.class)
  if (x$alternative=="greater") cat(" greater than ")
  if (x$alternative=="less") cat(" less than ")
  if (x$alternative=="two.sided") cat(" not equal to ")
  cat("the margins", "\n")

comparison <- rep(x$comp.names, each=length(x$resp))
endpoint <- rep(x$resp, times=length(x$comp.names))
margin <- estimate <- statistic <- p.value.raw <- p.value.adj <- NULL
for (i in 1:length(x$comp.names)) {
  margin <- c(margin, round(x$Margin[i,],digits))
  estimate <- c(estimate, round(x$estimate[i,],digits))
  statistic <- c(statistic, round(x$statistic[i,],digits))
  p.value.raw <- c(p.value.raw, round(x$p.val.raw[i,],digits))
  p.value.adj <- c(p.value.adj, round(x$p.val.adj[i,],digits))
}
out <- data.frame(comparison, endpoint, margin, estimate, statistic, p.value.raw, p.value.adj)
cat("", "\n")
print(out, digits=digits)
cat("", "\n")

invisible(x)


}
