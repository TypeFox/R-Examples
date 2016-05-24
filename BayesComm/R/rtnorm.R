rtnorm <-
function(n, mu, si, low, up){
	.Call( "rtnorm", n, mu, si, low, up, PACKAGE = "BayesComm" )
}