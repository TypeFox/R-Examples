
# Figure 6.5-1
xyplot(dbinom(0:18,18,0.5) ~ 0:18, type='h', 
	groups=0:18 %in% 0:12, lwd=8, lty=1)

plotDist('binom', params=c(18,0.5), groups=x %in% 0:12, kind='histogram') 

