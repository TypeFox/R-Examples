# tree2str.r #############################################################################################################
# FUNCTION:             	DESCRIPTION:
#  get.params				Prints the non-sorted parameter values by default.
#  .read.params 			Returns the parameter values from a copula tree. (Internal function)
#  tree2str					Prints the structure of HACs as string of class character.
#  .allocate.all			Returns and constructs the tree, which is printed by tree2str. (Internal function)
#  .one.with.theta	        Adds a leave to the string and calls .allocate.all with theta = TRUE. (Internal function)
#  .one.without.theta	    Adds a leave to the string and calls .allocate.all with theta = FALSE. (Internal function)
##########################################################################################################################

get.params = function(hac, sort.v = FALSE, ...){
	res = numeric(1)
    res = .read.params(hac$tree)
    if(sort.v == FALSE){res}else{sort(res, ...)}
}

#-------------------------------------------------------------------------------------------------------------------------------

.read.params = function(tree){
	if(length(tree)==1){tree = tree[[1]]}
	n = length(tree)
	s = sapply(tree, is.list)
		
		if(any(s)==TRUE){
			res = c(tree[[n]], unlist(sapply(X = tree[s], FUN = .read.params)))
		}else{
			res = tree[[n]]}
	res
}

#-------------------------------------------------------------------------------------------------------------------------------

tree2str = function(hac, theta = TRUE, digits = 2){
	res = character(1)
    res = .allocate.all(hac$tree, theta, digits)
    res
}

#-------------------------------------------------------------------------------------------------------------------------------

.allocate.all = function(tree, theta, digits){
n = length(tree); x = character(1)
	
		if(theta){
			for(i in 1:n){
				if(i == 1){
					x = paste("(", .one.with.theta(tree[[i]], digits = digits, theta = theta), sep = "")
			}else{
				if((i > 1) & (i < n)){
					x = paste(x, ".", .one.with.theta(tree[[i]], digits = digits, theta = theta), sep = "")
			}else{
				if(i == n){
					x = paste(x, ")_{", .one.with.theta(tree[[i]], digits = digits, theta = theta),"}", sep = "")
			}}}}
	}else{
		if(!theta){
			for(i in 1:(n-1)){
				if(i == 1){
					x = paste("(", .one.without.theta(tree[[i]], theta = theta), sep = "")
			}else{
				if((i > 1) & (i < n)){
					x = paste(x, ".", .one.without.theta(tree[[i]], theta = theta), sep = "")
			}}}}
					x = paste(x, ")", sep = "")
			}
	return(x)
}

#-------------------------------------------------------------------------------------------------------------------------------

.one.with.theta = function(element, digits, theta){
		if(class(element)=="character"){
			d = length(element)
			if(d ==1){
				element
		}else{
			if(d > 1){
				paste(element, sep = "", collapse = ".")
		}}
	}else{
		if(class(element)=="list"){
			.allocate.all(element, theta = theta, digits = digits)
	  }else{
			round(element, digits = digits)
	}}
}

#-------------------------------------------------------------------------------------------------------------------------------

.one.without.theta = function(element, theta){
		if(class(element)=="character"){
			d = length(element)
			if(d ==1){
				element
		}else{
			if(d > 1){
				paste(element, sep = "", collapse = ".")
		}}
	}else{
		if(class(element)=="list"){
			return(.allocate.all(element, theta = theta))
	}}
}
