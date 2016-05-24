setClassUnion("list.null", c("list","NULL"))
.Valid.splineMatrix <- function(object) {

	errors <- NULL

	tmp.test.null <- 
		all(!is.null(object@Knots) & 
			!is.null(object@Boundaries) & 
			!is.null(object@Content_Areas) & 
			!is.null(object@Grade_Progression) & 
			!is.null(object@Time) & 
			!is.null(object@Time_Lags) & 
			!is.null(object@Version))

	tmp.test.length <- length(object@Content_Areas[[1]]) == length(object@Grade_Progression[[1]]) & length(object@Grade_Progression[[1]]) == length(object@Time[[1]])

	if (!is.null(object@Content_Areas) & !is.null(object@Grade_Progression) & !is.null(object@Time)) {
		tmp.test.character <- is.character(object@Content_Areas[[1]]) & is.character(object@Grade_Progression[[1]]) & is.character(object@Time[[1]])
	}

	if (!is.null(object@Time)) {
		if (any(object@Time[[1]]=="BASELINE")) {
			tmp.test.time <- TRUE
		} else {
			tmp.test.time <- all(yearIncrement(head(object@Time[[1]], 1), c(0, cumsum(object@Time_Lags[[1]]))) == object@Time[[1]])
		}
	}

	if (!is.null(object@Time_Lags)) {
		tmp.test.numeric.lags <- is.numeric(object@Time_Lags[[1]])
	}
	
	if (!tmp.test.null) errors <- c(errors, "\tError in splineMatrix construction: Slots in splineMatrix are not all non-null.")
	if (!tmp.test.length) errors <- c(errors, "\tError in splineMatrix construction: Slots @Grade_Progression, @Content_Areas, and @Time are not all the same length.")
	if (!tmp.test.character) errors <- c(errors, "\tError in splineMatrix construction: Slots @Grade_Progression, @Content_Areas, and @Time are not all characters.")
	if (!tmp.test.time) errors <- c(errors, "\tError in splineMatrix construction: Slots @Time and @Time_Lags are not commensurate with one another.")
	if (!tmp.test.numeric.lags) errors <- c(errors, "\tError in splineMatrix construction: Slot @Time_Lags is not numeric.")

	if (is.null(errors)) TRUE else errors
}
setClass("splineMatrix", 
	contains='matrix', 
	representation(
		Knots='list.null', 
		Boundaries='list.null', 
		Content_Areas='list.null', 
		Grade_Progression='list.null', 
		Time='list.null', 
		Time_Lags='list.null', 
		Version="list.null"),
		validity=.Valid.splineMatrix)
setValidity("splineMatrix", .Valid.splineMatrix)
