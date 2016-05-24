find.base1 <- function(time, fixations, conc) {
	if(conc) {
		min(fixations[time < find.mu(time, fixations, conc)])
	} else max(fixations[time < find.mu(time, fixations, conc)])
}