# Automatic creation of dense two-way classification design matrix
# for usage of dense robust estimation with rq.fit.sfn (package quantreg).
# or slm.fit (package SparseM).
# The sum of the second factor is constrained to be zero. No general mean.
# Usually this is much easier created by:
# >D[,2] <- factor(D[,2])
# >ny <- length(levels(D[,2]))
# >D[,3] <- factor(D[,3])
# >ns <-  length(levels(D[,3]))
# >o <-append(D[,1],0)
# >X1 <- as.matrix.csr(model.matrix(~D[,2]-1))
# >X2 <- as.matrix.csr(model.matrix(~D[,3]-1))
# >X3 <- as.matrix.csr(t(c(rep(0,ny),rep(1,ns))))
# >MM <- rbind(cbind(X1,X2),X3)
# >slm.fit(MM,o)$coef
# this gives something similar:
# ddm <- as.matrix.csr(model.matrix(~ y + s -1, contrasts=list(s=("contr.sum"))))
# however, this procedure is quite memory demanding and might exceed storage
# capacity for large problems. Moreover, in order to get direct estimates for all
# coefficients, an additional row appended to the matrix, where the colums
# for the second factor are set to 1.
# This procedure here is much less memory comsuming.
# input: data frame with three columns: (observations, factor 1, factor 2)
# (phenlogy: observation day of phase, year, station)
# output: dense roworder matrix, matrix.csr format (see matrix.csr in package SparseM)
# and the sorted data frame D (data frame is being sorted first by f2 then by f1 )
# rows with NA's are removed by default
pheno.ddm <- function(D, na.omit=TRUE) {
	if(!is.data.frame(D) || length(D) != 3) 
		stop("phen.ddm: argument must be data frame with 3 fields. Exiting ...")

	# rows with NA's
	rows.na <- which(is.na(D)==TRUE)%%length(D[[1]])
	if(na.omit & length(rows.na) > 0) {
		D <- D[-rows.na,]
		# order first by factor 2 then by factor 1
		D <-  D[order(D[[3]],D[[2]]),]
		n <- length(D[[1]])	 	# number of observations
		f1 <- factor(D[[2]]) 		# factor 1: year
		n1 <- nlevels(f1) 		# number of levels factor 1 (phenlogy: years)
		f2 <- factor(D[[3]])		# factor 2: station
		n2 <- nlevels(f2)		# number of levels factor 2 (phenlogy: station)
	}
	else {
		# order first by factor 2 then by factor 1
		D <-  D[order(D[[3]],D[[2]]),]
		n <- length(D[[1]])	 	# number of observations
		f1 <- factor(D[[2]]) 		# factor 1: year
		n1 <- nlevels(f1) 		# number of levels factor 1 (phenlogy: years)
		f2 <- factor(D[[3]])		# factor 2: station
		n2 <- nlevels(f2)		# number of levels factor 2 (phenlogy: station)
	}
	
	# ra: Object of class numeric, a real array of nnz elements containing the 
	#	nn-zero elements of A, stored in row order. Thus, if i<j, all elements 
	#	of row i precede elements from row j. The order of elements within the 
	#	rows is immaterial. 
	# ja: Object of class integer, an integer array of nnz elements containing 
	#	the column indices of the elements stored in  ra. 
	# ia: Object of class integer, an integer array of n+1 elements containing 
	#	pointers to the beginning of each row in the arrays  ra  and  ja . 
	#	Thus  ia[i]  indicates the position in the arrays  ra  and  ja  where 
	#	the ith row begins. The last, (n+1)st, element of  ia  indicates 
	# 	where the n+1 row would start, if it existed. 
	# dimension: Object of class integer, dimension of the matrix
	
	# number of nn-zero elements in ra and ja.
	# each observation has an entry for year and station,
	# additionally the contraint that the sum of station effects is zero add n2 
	nnz <- n*2+n2
	ra=numeric(nnz)
	ja=integer(nnz)
	ia=integer(n+2)

	dimension=integer(2)
	ra[1:nnz] <- 1
	ia[1] <- as.integer(1) 
	# the observations
	for(i in 1:n) {
		ja[2*i-1] <- as.integer(f1[i])
		ja[2*i] <- n1 + as.integer(f2[i])
		ia[i+1] <- as.integer(ia[i] + 2)
	}
	# the sum of the station effects
	ia[n+2] <- as.integer(nnz+1)
	for(i in 1:n2) {
		ja[n*2+i] <- n1 + as.integer(i)
	}

	dimension <- as.integer(c(n+1,n1+n2))
	ddm <- new("matrix.csr",ra=ra,ja=ja,ia=ia,dimension=dimension)
	return(list(ddm=ddm,D=D,rows.na=rows.na))
}
