`summary.etc.rat` <-
function(object,digits=4,...) {


cat("", "\n")
cat("Alternative hypotheses: ratios to control within the margins", object$margin.lo, "and", object$margin.up,
    "\n")
cat("Method:", object$method,
    "\n")
cat("", "\n")

if (object$method=="var.equal") {
  cat("Degree of freedom:", "\n")
  print(object$degr.fr, digits=digits)
  cat("", "\n")
  cat("Critical value:", "\n")
  print(object$crit.value, digits=digits)
  cat("", "\n")
  if (is.numeric(object$conf.int)) {
    out <- cbind(object$estimate, object$test.stat, object$conf.int["lower",], object$conf.int["upper",],
                 object$p.value)
    colnames(out) <- c("estimate", "statistic", "lower", "upper", "p.value")
  } else {
    out <- cbind(object$estimate, object$test.stat, object$p.value)
    colnames(out) <- c("estimate", "statistic", "p.value")
    cat("                                      ", "\n")
    cat("   The mean in the denominator is not significantly different from zero. ", 
        "\n")
    cat("                                      ", "\n")
  }
}

if (object$method=="var.unequal") {
  cat("Degree of freedom:", "\n")
  print(rbind(up=object$degr.fr.up,down=object$degr.fr.do,SCI=object$degr.fr.ci), digits=digits)
  cat("", "\n")
  cat("Critical value:", "\n")
  print(rbind(up=object$crit.value.up,down=object$crit.value.do,SCI=object$crit.value.ci), digits=digits)
  cat("", "\n")
  if (is.numeric(object$conf.int)) {
    out <- cbind(object$estimate, object$test.stat.up, object$test.stat.do,
                 object$conf.int["lower",], object$conf.int["upper",], object$p.value)
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

