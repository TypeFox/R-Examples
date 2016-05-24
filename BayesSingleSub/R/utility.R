logMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogMeanExpLogs", as.numeric(v), N, package="BayesSingleSub")
}