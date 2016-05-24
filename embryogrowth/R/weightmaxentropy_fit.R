

.weightmaxentropy_fit <- function(par, pi, entropy.method, pentropy, vm) {

result <- NULL

vec.weight <- c(abs(par), 1)
names(vec.weight)[length(vec.weight)] <- vm

  nrow <- dim(pi)[1]
  ncol <- dim(pi)[2]
  
  # juste pour reprendre la structure
	pi2 <- pi

  for(i in 1:nrow)
	pi2[i,] <- pi[i, 1:ncol]*(vec.weight[i]/sum(vec.weight))
  
  piglobal <- NULL
  for(i in 1:ncol) piglobal <- c(piglobal, sum(pi2[,i]))
  
	L <- modifyList(list(y=piglobal), pentropy)
	E <- do.call(entropy.method, L) 
  
	return(as.numeric(E))
  
}
