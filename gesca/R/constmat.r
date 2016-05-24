constmat <- function (A0)
{
	#---------------------------------------
	# generate orthogonal projector matrix for equality constraints on A
	# Heungsun Hwang, Sunmee Kim
	# Last revised Aug 23, 2015
	#---------------------------------------
	
	sizea <- dim(A0)
	vect0 <- matrix(A0, nrow = sizea[1]*sizea[2], byrow=T)
	vect <- vect0[vect0 >= 1] 	# all free parameters and equality constraints in A across groups
	nzct <- which(vect0 >= 1)
	nzt <- length(nzct)			# number of free (constrained) parameters in T across groups
	nzcst <- which(vect0 >= 1 & vect0 != 99) # all constrained parameters in T across groups
	num_nzct <- length(nzcst)	# number of constrained parameters in T across groups
	
	if (num_nzct == 0) {
		PHT <- diag(1,nzt)		# PHT = Projector of H
		num_const <- 0
	} else {
		const_t <- vect0[vect0 >= 1 & vect0 != 99]	# constraint values
		num_const <- max(const_t)					# number of constraint sets
		PHT <- diag(1,nzt)
		for (j in 1:num_const) {
			cont <- which(vect == j)
			num_cont <- length(cont)
			pht <- matrix(0,nzt,1)
			for (i in 1:num_cont) {
				pht[cont[i]] <- 1/num_cont
			}
			for (s in 1:num_cont) {
				PHT[,cont[s]] <- pht
			}
		}
	}	
	
	output.constmat <- list(PHT = PHT, num_nzct = num_nzct, num_const = num_const)
	output.constmat
}