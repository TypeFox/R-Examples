multcor.test <-
function(x, n = NULL, Df = NULL,
	alternative = c("two.sided", "less", "greater"),
	adjust = "none")
{
    stopifnot(is.matrix(x))
    if (nrow(x) != ncol(x))
       stop("'x' must be a square matrix!")
    if (any(abs(x) > 1))
       stop("'x' contains values less than -1 or greater than 1!")
    p <- ncol(x)
    if (is.null(n) & is.null(Df))
       stop("one of 'n' or 'Df' must be specified!")
    if (!is.null(n)) Df <- n - 2

   se <- sqrt((1 - x^2) / Df)
   tval <- x / se

   alternative <- match.arg(alternative, c("two.sided", "less", 
      "greater"))
   if (alternative == "two.sided") {
      pval <- 2 * pt(abs(tval), Df, lower.tail = FALSE)
      } else if (alternative == "less") {
      pval <- pt(tval, Df, lower.tail = TRUE)
      } else if (alternative == "greater") {
      pval <- pt(tval, Df, lower.tail = FALSE)
    }
   if (adjust != "none") {
      low.p <- pval[lower.tri(pval)]
      pval.adj <- p.adjust(low.p, adjust)
      pval[lower.tri(pval)] <- pval.adj
      pval <- as.matrix(as.dist(pval))
   }
   diag(pval) <- NA
   pval. <- round(pval, 4)
   for(i in 1:p) {
      for(j in 1:p) {
      if (i < j) pval.[i, j] <- indicate.signif(pval.[i, j])
      }
   }

   out <- list(t.values = tval, p.values = pval,
      p.check = noquote(pval.),
      adjustment = adjust, df = Df,
	alternative = alternative,
      data.name = deparse(substitute(x)))
   class(out) <- "multcor.test"
   return(out)
}

# ------------------------------------------------------
# print method
print.multcor.test <- 
function (x, digits = 4L, quote = TRUE, ...) 
{
    cat("\n")
    cat(strwrap("Pairwise correlation t-test", prefix = "\t"), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")
    if (!is.null(x$df)) 
       cat("degrees of freedom:", x$df, "\n")
    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        alt.char <- switch(x$alternative, two.sided = "not equal to 0", 
           less = "less than 0", greater = "greater than 0")
        cat("the true correlation is", alt.char, "\n")
    }
    if (!is.null(x$p.check)) {
        cat("p-values (with", x$adjustment, "adjustment for multiple tests): \n")
        print(x$p.check, ...)
        cat("---",
        "\nSignif. codes: '***'0.001 '**'0.01 '*'0.05 '.'0.1 ' '1 \n")
    }
    cat("\n")
    invisible(x)
}

