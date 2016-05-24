`summary.etc.diff` <-
function(object,digits=4,...) {


cat("", "\n")
cat("Alternative hypotheses: differences to control within the margins", object$margin.lo, "and", object$margin.up,
    "\n")
cat("Method:", object$method,
    "\n")
cat("", "\n")
if (object$method!="non.par")
{
  cat("Degree of freedom:", "\n")
  print(object$degr.fr, digits=digits)
  cat("", "\n")
  cat("Critical value:", "\n")
  print(object$crit.value, digits=digits)
  cat("", "\n")
}
if (object$method=="Bofinger")
{
  cat("Correlation matrix:", "\n")
  print(object$corr.mat, digits=digits)
  cat("", "\n")
}

if (object$method!="non.par")
{
  out <- cbind(object$estimate, object$test.stat, object$conf.int["lower",], object$conf.int["upper",],
               object$p.value)
  colnames(out) <- c("estimate", "statistic", "lower", "upper", "p.value")
} else {
  out <- cbind(object$estimate, object$test.stat[,"test.stat.up"], object$test.stat[,"test.stat.do"], object$conf.int["lower",],
               object$conf.int["upper",], object$p.value)
  colnames(out) <- c("estimate", "statistic.up", "statistic.do", "lower", "upper", "p.value")
}


print(out, digits=digits)
cat("", "\n")


}

