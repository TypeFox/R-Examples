
# Figure 6.4-1
xyplot(dbinom(0:18,27,0.25) ~ 0:18, groups=0:18 %in% 7, 
	type='h', lwd=8, lty=1)
plotDist('binom', params=c(18,.5), groups= x %in% 7, kind='histogram') 

