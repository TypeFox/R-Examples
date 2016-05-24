rosenberg <- function(phy){	
	RosenbergP_AB <- function(a, b){
		d1 <- choose(a+b, a)
		d2 <- a+b-1
		n1 <- 2
		n2 <- 1
		(n1/d1)*(n2/d2)
	}
	mat <- polyBalance(phy)
	lab <- apply(mat, 1, function(x) RosenbergP_AB(x[1], x[2]))
	lab
}

