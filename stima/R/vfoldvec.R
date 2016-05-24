vfoldvec <-
function(n,k=10,seed=3){
	    #k equals number of sets
	   nset <- trunc(n/k)
		rest <- n - nset * k
		set.seed(seed)
		vec <- sample(1:k)
		if(nset != 1) {
				for(i in 1:(nset - 1)) {
				vec <- c(vec, sample(1:k))
				}
		}
		if(rest != 0) {
			vectot <- c(vec[1:floor(nset/2 * k)], sample(1:rest), vec[floor(
			(nset/2 * k) + 1):(nset * k)])
			}
			else {vectot <- vec}
				return(vectot)}
