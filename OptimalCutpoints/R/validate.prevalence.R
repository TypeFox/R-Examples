validate.prevalence <-
function(prev) {
	if (prev < 0 || prev > 1) {
		stop("You have entered an invalid value for prevalence. \n Prevalence must be a value higher than 0 and lower than 1.", call. = FALSE)
	}
	if (prev == 0) {
		stop("You have entered an invalid value for prevalence. \n No subject in the population has the disease. Please check this value and \n introduce another valid value. \n Prevalence must be a value higher than 0 and lower than 1.", call. = FALSE)
	}
	if(prev == 1) {
		stop("You have entered an invalid value for prevalence. \n All subjects in the population have the disease. Please check this value and \n introduce another valid value. \n Prevalence must be a value higher than 0 and lower than 1.", call. = FALSE)
  	}
}
