# Publicly facing quick functions for calculating jbar, jstar, and pstar 

##########
# jbar   #
##########

jbar = function(x){
	mean(1 - vegdist(t(x), method = "jaccard"))
}



##########
# jstar  #
##########

jstar = function(x, n = NULL){
	UseMethod("jstar")
}

jstar.matrix = function(x, n = NULL){
	jstar(as.data.frame(x))
}

jstar.data.frame = function(x, n = NULL){
	
	if(!missing(n)){
		if(n != ncol(x)) warning("n does not match ncol(x)")
	}
	jstar(rowMeans(x), ncol(x))
}

jstar.numeric = function(x, n){
	
	if(any(x > 1) | any(x < 0)){
		stop("occupancy rates must be between 0 and 1")
	}
	
	# "round" is necessary to prevent a floating point error.
	`/`(
		sum(choose(round(x * n), 2)),	
		sum(choose(n, 2) - choose(round((1 - x) * n), 2))
	)
}

#########
# pstar #
#########

pstar = function(x, n = NULL){
	UseMethod("pstar")
}

pstar.matrix = function(x, n = NULL){
	pstar(as.data.frame(x))
}

pstar.data.frame = function(x, n = NULL){
	
	if(!missing(n)){
		if(n != ncol(x)) warning("n does not match ncol(x)")
	}
	pstar(rowMeans(x), ncol(x))
}


pstar.numeric = function(x, n){
	j = jstar(x, n)
	
	`/`(
		j * (2 * n - 1) + 1,
		(j+1) * n
	)
}

