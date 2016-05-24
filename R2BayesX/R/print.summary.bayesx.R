print.summary.bayesx <- function(x, digits = max(3, getOption("digits") - 3),
  signif.stars = getOption("show.signif.stars"), ...)
{
  n <- length(x)
  ncheck <- n > 1L
  digits <- attr(x, "digits")
  nx <- names(x)
  if(is.null(nx))
    nx <- 1L:n
  for(i in 1L:n) {
    if(ncheck)
      cat("###", nx[i], "\n")
    .print_summary_bayesx(x[[i]], digits = digits, signif.stars = signif.stars, ...)
  }
  if(ncheck) {
    cat("###\n")
    cat("Object consists of", n, "models\n")
  }

  return(invisible(NULL))
}

