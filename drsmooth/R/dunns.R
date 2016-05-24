#' @rdname dunns
#' @title Dunn's Test
#' @param targetcolumn  Character string, name of response column to be tested.
#' @param alternative  The direction(s) of the alternative hypothesis.
#' @param alpha  Significance level (numeric) to be used.
#' @param control  Not used in this test.
#' @param tot.obs  Not used in this test.
#' @param label Label indicating alternative.
#' @param data  Input dataframe.
#' @keywords internal

dunns <- function (targetcolumn, alternative, alpha, control, tot.obs, label, data) {
  if (alternative == "greater") {cont <- "one-tailed"} else {
    if (alternative == "two.sided") {cont <- "two-tailed"}
  }
	kruskal <- pgirmess::kruskalmc(data[,targetcolumn], data$dose_fac, probs=alpha, cont = cont)
  title <- paste("Dunn's", label, "Test", sep = " ")
  
  # R CMD check resists internal ::: call, but :: fails
  #output <- drsmooth:::dunns.format(kruskal, title)
  f <- get("dunns.format", envir = environment(drsmooth))
  output <- f(kruskal, title)
  
	return(output)
}