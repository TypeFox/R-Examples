# Pr[ 3 males ]
dbinom(3,5,0.2)
# or do it by hand
choose(5,3) * (0.2)^3 * (1-0.2)^(5-3)

# Pr[ 6 left-handed flowers ]
dbinom(6, 27, 0.25)
choose(27,6)

# Table 7.1-1
cbind( x=0:27, prob=dbinom(0:27, 27, 0.25) ) 
# Figure 7.1-1
xyplot(dbinom(0:27, 27, 0.25) ~ 0:27, type='h', lwd=4)
histochart(dbinom(0:27, 27, 0.25) ~ 0:27)
# Figure 7.1-2
xyplot(dbinom(0:10, 10, 0.25) ~ (0:10)/10, type='h', lwd=6, 
	xlab='proportion of successes',
	ylab='probability'
	)
xyplot(dbinom(0:100, 100, 0.25) ~ (0:100)/100, type='h', lwd=4, 
	xlab='proportion of successes',
	ylab='probability'
	)
histochart(dbinom(0:10, 10, 0.25) ~ (0:10)/10, 
	xlab='proportion of successes',
	ylab='probability',
	scales=list(x=list(draw=FALSE))
	)
histochart(dbinom(0:100, 100, 0.25) ~ (0:100)/100, 
	xlab='proportion of successes',
	ylab='probability',
	scales=list(x=list(draw=FALSE))
	)

