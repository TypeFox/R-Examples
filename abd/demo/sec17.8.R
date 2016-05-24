
# Figure 17.8-3 but using LOESS rather than splines
xyplot(length ~ age, ShrinkingSeals, cex=.5, alpha=.8, col.line='black',
	type=c('p','smooth'), span=.1, degree=2, eval=500
	)
