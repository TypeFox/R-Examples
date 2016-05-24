summary.multeq.rat <-
function(object,digits=4,...) {


cat("", "\n")
cat("Alternative hypotheses: ratios")
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

if (object$var.equal==TRUE) {
  cat("Degree of freedom:", "\n")
  print(object$degr.fr, digits=digits)
  cat("", "\n")
  if (is.numeric(object$lower) & is.numeric(object$upper)) {
    out <- cbind(object$estimate, object$test.stat, object$lower, object$upper, object$p.value)
    colnames(out) <- c("estimate", "statistic", "lower", "upper", "p.value")
  }
  else {
    out <- cbind(object$estimate, object$test.stat, object$p.value)
    colnames(out) <- c("estimate", "statistic", "p.value")
    cat("                                      ", "\n")
    cat("   The mean in the denominator is not significantly different from zero. ", 
        "\n")
    cat("                                      ", "\n")
  }
}

if (object$var.equal==FALSE) {
  cat("Degree of freedom:", "\n")
  print(rbind(up=object$degr.fr.up,down=object$degr.fr.do,SCI=object$degr.fr.ci), digits=digits)
  cat("", "\n")
  if (is.numeric(object$lower) & is.numeric(object$upper)) {
    out <- cbind(object$estimate, object$test.stat.up, object$test.stat.do,
                 object$lower, object$upper, object$p.value)
    colnames(out) <- c("estimate", "statistic(up)", "statistic(down)", "lower", "upper", "p.value")
  } else {
    out <- cbind(object$estimate, object$test.stat.up, object$test.stat.do, object$p.value)
    colnames(out) <- c("estimate", "statistic.up", "statistic.do", "p.value")
    cat("                                      ", "\n")
    cat("   The mean in the denominator is not significantly different from zero. ", 
        "\n")
    cat("                                      ", "\n")
  }
}

print(out, digits=digits)
cat("", "\n")


}

