# aggregate.r ############################################################################################################
# FUNCTION:         DESCRIPTION:
#  .union         	Aggregates the whole tree of a HAC. (Internal function)
#  .compare.one     Compares two paramters of subcopulae for being larger or smaller than the prespecified epsilon and aggregates the subcopulae if the condition is TRUE. (Internal function)
#  .method					Computates the mean, min and max if abs(theta_1 - theta_2) < epsilon. (Internal function)
#  aggregate.hac    Aggregates an object of the class 'hac'.             
##########################################################################################################################

.union = function(tree, epsilon, method, ...){
	n = length(tree)
	result = list()
	
	for(i in 1:n){
	if(class(tree[[i]])=="list"){
			x1 = list()
			x2 = list()
			
		if(i == 1){
			for(j in 1:(n-i)){x2[[j]]=tree[[j+i]]}
			tree.without.i = c(x2)
		}else{
		if(i > 1){
			for(j in 1:(i-1)){x1[[j]]=tree[[j]]}
			for(j in 1:(n-i)){x2[[j]]=tree[[j+i]]}
			tree.without.i = c(x1, x2)
		}}
		
		tree = .compare.one(tree.one = tree.without.i, tree.two = tree[[i]], epsilon = epsilon, method = method, ...)
		n = length(tree)
	}}
	tree
}
 
#--------------------------------------------------------------------------------------------

.compare.one = function(tree.one, tree.two, epsilon, method, ...){
			n1 = length(tree.one)
			n2 = length(tree.two)
			
		if(abs(tree.one[[n1]]-tree.two[[n2]]) < epsilon){
			x1 = list()			
			x2 = list()
			
			for(j in 1:(n1-1)){x1[[j]]=tree.one[[j]]}
			for(j in 1:(n2-1)){x2[[j]]=tree.two[[j]]}
			
			tree.new = c(x1, x2, .method(c(tree.one[n1][[1]], tree.two[n2][[1]]), method = method, ...))
			.union(tree = tree.new, epsilon = epsilon, method = method, ...)
	}else{
			c(list(.union(tree.two, epsilon = epsilon, method = method, ...)), tree.one)
	}
}

#--------------------------------------------------------------------------------------------

.method = function(x, method, ...){
	if(method == "mean"){mean(x, ...)}
	else
	if(method == "min"){min(x, ...)}
	else
	if(method == "max"){max(x, ...)}
}

#--------------------------------------------------------------------------------------------

aggregate.hac = function(x, epsilon = 0, method = "mean", ...){
 	tree = .union(x$tree, epsilon = epsilon, method = method, ...)
 	hac(type = x$type, tree = tree)
}
