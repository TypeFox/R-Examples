distribution_estimate <-
function(v, successes=NULL, num_quantiles=101, observed=FALSE) {
	p_to_examine = seq(0,1,length.out=num_quantiles)
	if (observed) {
		if (!is.null(successes)) {
			total = v
			observed_p = successes/total
			quantiles = sapply(p_to_examine, function(p) {as.numeric(p > observed_p)})
		} else {
			quantiles = as.vector(quantile(v,p_to_examine,type=2))
		}
	} else {
		if (!is.null(successes)) {
			total = v
			quantiles = qbeta(p_to_examine,successes+1,total-successes+1)
		} else {
			bootstrapped_means = boot(v, function(y,i) {mean(y[i])}, 1000)
			quantiles = as.vector(quantile(bootstrapped_means$t,p_to_examine,type=2))
		}
	}
	x = rep(quantiles,each=2)
	mids = (quantiles[2:length(quantiles)] + quantiles[1:(length(quantiles)-1)])/2.0
	widths = quantiles[2:length(quantiles)] - quantiles[1:(length(quantiles)-1)]
	heights = 1/widths
	y = c(0,rep(heights,each=2),0)

	x=x[is.finite(y)]
	y=y[is.finite(y)]

	mids=mids[is.finite(heights)]
	widths=widths[is.finite(heights)]
	lefts = quantiles[1:(num_quantiles-1)]
	lefts=lefts[is.finite(heights)]
	rights=quantiles[2:num_quantiles]
	rights=rights[is.finite(heights)]
	probabilities = cumsum(rep(p_to_examine[2],num_quantiles-1))
	probabilities = probabilities[is.finite(heights)]
	probabilities = probabilities - c(0,probabilities[1:(length(probabilities)-1)])
	heights=heights[is.finite(heights)]
	heights = heights * (probabilities/p_to_examine[2])

	return(list(quantiles=quantiles,x=x,y=y,lefts=lefts,mids=mids,rights=rights,heights=heights,widths=widths,probabilities=probabilities))
}
