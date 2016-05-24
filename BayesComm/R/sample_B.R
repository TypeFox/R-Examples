sample_B <-
function (z, X) {
	.Call( "sample_B", z, X, PACKAGE = "BayesComm" )
}