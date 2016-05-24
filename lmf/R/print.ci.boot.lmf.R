print.ci.boot.lmf <-
function(x,
                              ...)
{
  cat("\nESTIMATING FLUCTUATING SELECTION IN AGE-STRUCTURED POPULATIONS\n",
      sep = "")
  cat(sprintf(ngettext(x$nboot, "%s Bootstrap replicate generated\n",
                       "%s Bootstrap replicates generated\n"), paste(x$nboot)), sep = "")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(sprintf("%s BOOTSTRAP CONFIDENCE INTERVALS:\n\n",
              paste((1 - x$clevel) * 100, '%', sep = "")), sep = "")
  if(!is.null(x$l))
  {
    cat("Projection matrix (l):\n")
    print.default(x$l, print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Lambda:\n", x$luv[1], "\n\n", sep = " ")
    nage <- dim(x$l)[1]
    cat("Stable age dist.(u) and reprod. values(v):\n")
    print.data.frame(data.frame(age = x$uage, u = x$luv[2 : (1 + x$nage)],
                                v = x$luv[(2 + x$nage) : (1 + x$nage + x$nage)]),
                     print.gap = 2, quote = FALSE, row.names = FALSE)
    cat("\n")
  }
  if(x$what == "alpha" | x$what == "all")
  {
    cat("Variance components:\n")
    cat(" Environmental\n", paste(" ", x$sigma2.e,
                                  sep = ""), "\n")
    cat(" Demographic\n")
    print.data.frame(data.frame(age = c(x$uage, "(total)"),
                                confint = c(x$sigma2.dj, x$sigma2.d)), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
    cat("Temporal mean alpha (a(M)):\n")
    print.default(x$aM, print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Temporal covariance matrix (M):\n")
    print.default(x$M, print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Temporal alpha estimates assuming no fluct. selection (a|M = 0):\n")
    print.default(x$anf, print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Covariance matrix assuming no fluct. selection (A):\n")
    print.default(x$Anf, print.gap = 2, quote = FALSE)
    cat("\n")
  }
  cat("-End-\n")
}
