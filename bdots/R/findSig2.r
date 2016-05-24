find.sig2 <- function(time, fixations, conc) {
	mu <- find.mu(time, fixations, conc)
	fixations <- fixations - find.base2(time, fixations, conc)
	fixations <- fixations[time >= mu]
	time <- time[time >= mu]
	total.fix <- sum(fixations)
	
	sigma <- time[which.min(abs((pnorm(1) - pnorm(-1)) * total.fix - cumsum(fixations)))] - mu
	sigma
}