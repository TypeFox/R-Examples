
alrt <- function(x1, x2, boundary = FALSE) {
  jll1 <- logLik(x1)
  jll2 <- logLik(x2)
  df1 <- attr(jll1, "df")
  df2 <- attr(jll2, "df")
  jll.diff <- abs(c(jll1) - c(jll2))
  df.diff <- abs(df1 - df2)
  p.value <- 1 - pchisq(2 * jll.diff, df = df.diff)
  if (boundary) p.value <- p.value / 2
  results <- list(out.tab = data.frame(model = c(1,2),
                                       jll = c(jll1, jll2),
                                       df = c(df1, df2)),
                  jll.diff = jll.diff,
                  df.diff = df.diff,
                  p = p.value)
  cat("\nLL of model 1: ", jll1, " df: ", df1, 
      "\nLL of model 2: ", jll2, " df: ", df2, 
      "\nDifference: ", jll.diff, " df: ", df.diff, 
      "\np-value against H_0: no difference between models ", 
      p.value, "\n")
  return(invisible(results))
}
