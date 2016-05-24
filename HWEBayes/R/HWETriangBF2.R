HWETriangBF2 <-
function(nvec){
	    if (length(nvec) != 3) stop("HWETriangBF2: Dimension of nvec not equal to 3\n")
	     bvecH1 <- c(1,1,1)
	     HWETriangBF2 <- TriangNormHWE(nvec)/DirichNormSat(nvec,bvecH1)
	     HWETriangBF2
}

