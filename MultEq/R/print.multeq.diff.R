print.multeq.diff <-
function(x,digits=4,...) {


cat("", "\n")
cat("Alternative hypotheses: differences")
if (is.numeric(x$margin.lo)) cat(" greater than", x$margin.lo)
if (is.numeric(x$margin.lo) & is.numeric(x$margin.up)) cat("\n", "and")
if (is.numeric(x$margin.up)) cat (" smaller than", x$margin.up)
cat("", "\n")
cat("Comparison:", x$comp.name,
    "\n")
cat("Method:", x$method,
    "\n")
cat("Equal variances:", x$var.equal,
    "\n")
out <- cbind(x$estimate, x$lower, x$upper, x$p.value)
colnames(out) <- c("estimate", "lower", "upper", "p.value")
cat("", "\n")
print(out, digits=digits)
cat("", "\n")
invisible(x)


}

