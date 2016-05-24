HWEDirichBF2 <-
function(nvec,bvec0,bvec1){
	      if (length(bvec0) != 2) stop("HWEDirichBF2: Dimension of bvec0 not equal to 2\n")
	      if (length(bvec1) != 3) stop("HWEDirichBF2: Dimension of bvec1 not equal to 3\n")
	      if (length(nvec) != 3) stop("HWEDirichBF2: Dimension of nvec not equal to 3\n")
	      w <- sum(bvec0); v <- sum(bvec1); n <- sum(nvec)
	      HWEDirichBF2 <- exp( nvec[2]*log(2) + lgamma(w) + lgamma(2*nvec[1]+nvec[2]+bvec0[1]) + lgamma(bvec1[1]) + lgamma(bvec1[2]) + lgamma(bvec1[3]) + lgamma(2*nvec[3]+nvec[2]+bvec0[2]) + lgamma(n+v) - lgamma(bvec0[1]) - lgamma(bvec0[2]) - lgamma(2*n+w) - lgamma(v) - lgamma(nvec[1]+bvec1[1]) - lgamma(nvec[2]+bvec1[2]) - lgamma(nvec[3]+bvec1[3]))
	     HWEDirichBF2 
}

