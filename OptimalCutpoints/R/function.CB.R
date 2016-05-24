function.CB <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc) {
	direction <- match.arg(direction)
	# The slope of the ROC curve at the optimal cutpoint is computed:
	S <- ((1-pop.prev)/pop.prev)*control$costs.ratio

	x <- (1 - measures.acc$Sp[,1])
	y <- measures.acc$Se[,1]

	rad <- (x^2 + y^2)^0.5
	theta <- atan2(y, x)

	theta.S <- atan(S) 
	theta.new <- theta - theta.S
	x.new <- rad * cos(theta.new)
	y.new <- rad * sin(theta.new)
 
	cCB <- measures.acc$cutoffs[which(round(y.new,10) == round(max(y.new, na.rm = TRUE), 10))]
				
	optimal.cutoff <- obtain.optimal.measures(cCB, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, optimal.criterion = S) 
	res
}
