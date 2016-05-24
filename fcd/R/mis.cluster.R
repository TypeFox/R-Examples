mis.cluster <- function(x, K, n){

	comms.list = list()
	true.list = list()
	for(j in 1: K){
		comms.list[[j]] = which(x == j)
		true.list[[j]] = c(((j-1)*n+1): (j*n))
	}
	
	list.perms = permn(1:K)
	N = length(list.perms)
	count = rep(0, N)
	for(i in 1: N){
		temp.perms = list.perms[[i]]
		for(j in 1: K){
			count[i] = count[i] + length(setdiff(true.list[[temp.perms[j]]], comms.list[[j]]))
		}
	}
	
	final = which.min(count)
	min.count = min(count)
	
	return(list(perms = list.perms[[final]], mis.cluster = min.count))
}