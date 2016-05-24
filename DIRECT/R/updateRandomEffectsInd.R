updateRandomEffectsInd <-
function (i, sdRET.curr, ts.cluster, cluster.mean, n.time, n.rep, upper, sdProp, DATA.FORMAT=0)
{
	tmp = rnorm (1, mean = sdRET.curr[i], sd = sdProp)
	sdRET.new = sdRET.curr
	sdRET.new[i] = reflect (tmp, lower=0, upper=upper)
	
	par.new = c (cluster.mean, sdRET.new)
	par.curr = c (cluster.mean, sdRET.curr)
	target.new = computeLogLik (par.new, ts.cluster, n.time, n.rep, DATA.FORMAT)
	target.curr = computeLogLik (par.curr, ts.cluster, n.time, n.rep, DATA.FORMAT)

	accept.prob.log = min (0, target.new - target.curr)
	if (log (runif(1)) < accept.prob.log)
		result = sdRET.new
	else
		result = sdRET.curr
	
	return (result)
}

