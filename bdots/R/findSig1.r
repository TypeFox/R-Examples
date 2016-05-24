find.sig1 <- function(time, fixations, conc) {
	mu <- find.mu(time, fixations, conc)
	fixations <- fixations - find.base1(time, fixations, conc)
	fixations <- rev(fixations[time <= mu])
	time <- rev(time[time <= mu])
	total.fix <- sum(fixations)
	
	sigma <- mu - time[which.min(abs((pnorm(1) - pnorm(-1)) * total.fix - cumsum(fixations)))]
	sigma
}