bdots.write.csv <- function(part2.list, file, ...) {
	group.1 <- as.character(part2.list$groups[1])
	group.2 <- as.character(part2.list$groups[2])

	out <- cbind(part2.list$curve.ci1, part2.list$curve.ci2, part2.list$sig)
	out <- out[,-5]
	colnames(out)[2:8] <- c(paste(group.1, "- Lower CI"), paste(group.1, "- Estimate"),
													paste(group.1, "- Upper CI"), paste(group.2, "- Lower CI"),
													paste(group.2, "- Estimate"), paste(group.2, "- Upper CI"),
													"Significance")
													
	write.csv(out, file, ...)
}