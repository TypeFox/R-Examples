drop.response <- function(formula, data) {
  tt <- terms(formula, data=data)
	if (length(attr(tt, "term.labels"))!=0) {
		formula <- reformulate(attr(tt, "term.labels"),
                           intercept = attr(tt, "intercept"))
	} else formula <- ~1
	return(formula)
}
