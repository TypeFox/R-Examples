find.mu <- function(time, fixations, conc) {
	if(conc) {
		time[which.max(fixations)]
	} else time[which.min(fixations)]
}