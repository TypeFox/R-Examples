"rbiased" <-
function(tip.number=20, p=0.3){
	
	tree <- matrix(NA, nrow = tip.number-1, ncol = 2)
	
	tree[tip.number-1,] <- c(-1,-2) 
	current.tip <- -3
	
	tip <- c(-1, -2)
	tab.prob <- c(p, 1-p)
	
	for ( node in (tip.number-2):1 ) {
		id <- sample(tip, 1, prob=tab.prob)
		tree[which(tree==id)] <- node
		tree[node,] <- c(id, current.tip)
		
		aux=sample(c(1,2),1)
		if (aux==1) {
			tab.prob <- c( tab.prob[tip!=id], tab.prob[tip==id]*c(p, 1-p))
		} else {
			tab.prob <- c( tab.prob[tip!=id], tab.prob[tip==id]*c(1-p, p))
		}
		
		tip <- c( tip[tip != id], id,  current.tip)
		current.tip <- current.tip-1
	}
	
	res <- treeshape(tree)
	res
}

