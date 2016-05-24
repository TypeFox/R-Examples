computeClusterParHastingsRatio <-
function (cluster.par.new, cluster.par.curr, gene.data, n.time, n.rep, DATA.FORMAT=0)
{
	loglik.new = computeLogLik (cluster.par.new, gene.data, n.time, n.rep, DATA.FORMAT)
	loglik.curr = computeLogLik (cluster.par.curr, gene.data, n.time, n.rep, DATA.FORMAT)
	
	logHRatio = loglik.new - loglik.curr
	
	return (logHRatio)
}

