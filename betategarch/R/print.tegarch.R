print.tegarch <-
function(x, ...)
{
vcovmat <- vcov(x)
out1 <- rbind(x$par, sqrt(diag(vcovmat)))
rownames(out1) <- c("Estimate:", "Std. Error:")
y.n <- length(x$y)
bic <- (-2*x$objective + length(x$par)*log(y.n))/y.n
out2 <- rbind(x$objective, bic)
rownames(out2) <- c("Log-likelihood:", "BIC:")
colnames(out2) <- ""

cat("Date:", x$date, "\n")
cat("Message (nlminb):", x$message, "\n")
if(!is.null(x$NOTE)){
  cat("NOTE:", x$NOTE, "\n")
}
cat("\n")
cat("Coefficients:\n")
print(out1)
print(out2)
}
