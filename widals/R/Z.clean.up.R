Z.clean.up <-
function(Z) {
	Z[ Z == Inf | Z == -Inf ] <- NA
	Z[ is.na(Z) | is.nan(Z) ] <- mean(Z, na.rm=TRUE)
	return(Z)
}
