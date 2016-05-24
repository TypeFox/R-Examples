rwish <-
function (S, df) {
	.Call( "rwish", S, df, PACKAGE = "BayesComm" )
}