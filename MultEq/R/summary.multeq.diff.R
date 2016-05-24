summary.multeq.diff <-
function(object,digits=4,...) {


cat("", "\n")
cat("Alternative hypotheses: differences")
if (is.numeric(object$margin.lo)) cat(" greater than", object$margin.lo)
if (is.numeric(object$margin.lo) & is.numeric(object$margin.up)) cat("\n", "and")
if (is.numeric(object$margin.up)) cat (" smaller than", object$margin.up)
cat("", "\n")
cat("Comparison:", object$comp.name,
    "\n")
cat("Method:", object$method,
    "\n")
cat("Equal variances:", object$var.equal,
    "\n")
cat("", "\n")
cat("Degree of freedom:", "\n")
print(object$degr.fr, digits=digits)
cat("", "\n")
out <- cbind(object$estimate, object$test.stat, object$lower, object$upper, object$p.value)
colnames(out) <- c("estimate", "statistic", "lower", "upper", "p.value")
print(out, digits=digits)
cat("", "\n")


}

