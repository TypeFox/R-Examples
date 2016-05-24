find.ht <- function(time, fixations, conc) {
	if(conc) {
		max(fixations)
	} else min(fixations)
}