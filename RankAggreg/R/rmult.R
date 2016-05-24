`rmult` <-
function(n=1, sspace, prob){
	prob <- prob/sum(prob)
	r <- runif(1)
	sspace[length(sspace)-sum(r < cumsum(prob))+1]
}
