rowSumsPrim <- function(x){
	.Call("R_rowSums", x
	,PACKAGE="gRbase"
	)
}

colSumsPrim <- function(x){
	.Call("R_colSums", x
	,PACKAGE="gRbase"
	)
}

colwiseProd <- function(vv, mm){
	.Call("R_colwiseProd", vv, mm
	,PACKAGE="gRbase"
	)
}
