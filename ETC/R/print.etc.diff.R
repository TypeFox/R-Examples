`print.etc.diff` <-
function(x,digits=4,...) {


cat("", "\n")
cat("Alternative hypotheses: differences to control within the margins", x$margin.lo, "and", x$margin.up,
    "\n")
cat("Method:", x$method,
    "\n")
out <- cbind(x$estimate, x$conf.int["lower",], x$conf.int["upper",], x$p.value)
colnames(out) <- c("estimate", "lower", "upper", "p.value")
cat("", "\n")
print(out, digits=digits)
cat("", "\n")
invisible(x)


}

