runif.faure <- function(n,dimension){

	pr <- c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,
			61,67,71,73,79,83,89,97,101, 103,107,109,
			113,127,131,137,139,149,151,157,163,167,173,
			179,181,191,193,197,199,211,223,227,229,233,
			239,241,251,257,263,269,271,277,281,283,293,
			307,311,313,317,331,337,347,349,353,359,367,
			373,379,383,389,397,401,409,419,421,431,433,
			439,443,449,457,461,463,467,479,487,491,499,
			503,509,521,523,541,547,557,563,569,571,577,
			587,593,599,601,607,613,617,619,631,641,643,
			647,653,659,661,673,677,683,691,701,709,719,
			727,733,739,743,751,757,761,769,773,787,797,
			809,811,821,823,827,829,839,853,857,859,863,
			877,881,883,887,907,911,919,929,937,941,947,
			953,967,971,977,983,991,997)

	# base of number system
	r <- pr[pr>=dimension][1]
	# compute coefficients
	# a priori upper bound for the length of the sequence
	m <- ceiling(log(n)/log(r))
	C <- matrix(0,m,m)
	C[1,] <- 1

	for (j in 2:m){
		C[j,j] <- 1      
		if ((j-1)>=2){
			for (i in 2:(j-1)){
				C[i,j] <- C[i,(j-1)] + C[(i-1),(j-1)]
			}
		}
	}
	C <- C%%r

	# to do the radical inverse
	dg <- r^(-(seq(1,m,1)))

	F  <- matrix(0,n,dimension)
	for (i in 1:n){
		# base r decomposition
		nr <- matrix(0,m,1)
		number <- i
		for (k in 1:m){
			quo <- floor(number/r)
			res <- number - quo*r
			number <- quo
			# i = nr(m)*r^(m-1) + ... + nr(2)*r + nr(1)
			nr[k] <- res 
		}
		F[i,1] <- sum(dg*nr)
		for (j in 2:dimension){
			nr <- (C%*%nr)%% r
			F[i,j] <- sum(dg*nr)
		}
	}
	# outputs
	return(out=list(design=F,n=n,dimension=dimension))
}