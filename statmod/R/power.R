power.fisher.test <- function(p1,p2,n1,n2,alpha=0.05,nsim=100,alternative="two.sided") {
#	Calculation of power for Fisher's exact test for
#	comparing two proportions
#	Gordon smyth
#	3 June 2003.  Revised 19 Nov 2011.

	y1 <- rbinom(nsim,size=n1,prob=p1)
	y2 <- rbinom(nsim,size=n2,prob=p2)
	y <- cbind(y1,n1-y1,y2,n2-y2)
	p.value <- rep(0,nsim)
	for (i in 1:nsim) p.value[i] <- fisher.test(matrix(y[i,],2,2),alternative=alternative)$p.value
	mean(p.value < alpha)
}
