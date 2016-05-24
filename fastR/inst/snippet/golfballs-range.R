stat <- function(x) { diff(range(x)) }  
plot <- statTally(golfballs,rgolfballs,stat,
	xlab="test statistic (range)")
