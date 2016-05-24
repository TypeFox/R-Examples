find.cross <- function(time, fixations) {
	med <- (max(fixations) - min(fixations)) * .5
	med.time <- time[which.min(abs(med - fixations))]
	med.time
}