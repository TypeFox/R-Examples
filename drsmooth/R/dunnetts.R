#' @name dunnetts
#' @title Dunnett's
#' @param targetcolumn  Character string, name of response column to be tested.
#' @param alternative Character string(s) specifying the direction of the alternative hypothesis.
#' @param alpha  Significance level (numeric) to be used.
#' @param control  Not relevant for this function
#' @param tot.obs  Not relevant for this function
#' @param label  Label of the alternative direction for output.
#' @param data  Input dataframe.
#' @keywords internal

dunnetts <- function (targetcolumn, alternative, alpha, control, tot.obs, label, data) {

	dose_fac <- data$dose_fac
    
	data.anova <- stats::aov(data[,targetcolumn] ~ dose_fac, data=data)
	dunnetts.glht <- multcomp::glht(data.anova, linfct=multcomp::mcp(dose_fac = "Dunnett"), alternative=alternative)
    
	summary <- summary(dunnetts.glht)
  levels <- levels(data$dose_fac)
  title <- paste("Dunnett's", label, "Multiple Comparisons Test", sep = " ")
  
	f <- get("dunnetts.format", envir = environment(drsmooth))
  output <- f(summary, title, levels)
  return(output)
}
