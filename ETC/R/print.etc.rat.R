`print.etc.rat` <-
function(x,digits=4,...) {


cat("", "\n")
cat("Alternative hypotheses: ratios to control within the margins", x$margin.lo, "and", x$margin.up,
    "\n")
cat("Method:", x$method,
    "\n")
if (is.numeric(x$conf.int)) {
  out <- cbind(x$estimate, x$conf.int["lower",], x$conf.int["upper",], x$p.value)
  colnames(out) <- c("estimate", "lower", "upper", "p.value")
}
else {
  out <- cbind(x$estimate, x$p.value)
  colnames(out) <- c("estimate", "p.value")
  cat("                                      ", "\n")
  cat("   The mean in the denominator is not significantly different from zero. ", 
      "\n")
  cat("                                      ", "\n")
}
cat("", "\n")
print(out, digits=digits)
cat("", "\n")
invisible(x)


}

