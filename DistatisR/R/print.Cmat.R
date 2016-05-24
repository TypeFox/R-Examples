print.Cmat <-
function (x,...) {

	res.Cmat <- x
	if (!inherits(res.Cmat, "Cmat")) stop ("no convenient data")
	cat("**Results related to the C matrix**\n")
	cat("*The results are available in the following objects:\n\n")
		
	if(res.Cmat$compact){
		#list(alpha=alpha) 
		res <- array("", c(1, 2), list(1:1, c("name", "description")))
		res[1,] <- c("$alpha","alpha weights")
	}else{
		#list(C = C, eigVector = eigC$vector, eigValues = eigC$values, tau = eigC$tau, G = eigC$G, alpha=alpha)		
		res <- array("", c(6, 2), list(1:6, c("name", "description")))
		res[1,] <- c("$C","The C matrix")
		res[2,] <- c("$eigVector","Eigenvectors")
		res[3,] <- c("$eigValues","Eigenvalues")
		res[4,] <- c("$tau", "Explained variance for components")
		res[5,] <- c("$G","the Rv matrix")
		res[6,] <- c("$alpha","alpha weights")
	}

	print(res)
}
