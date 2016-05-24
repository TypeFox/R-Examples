# constructor.r ##########################################################################################################
# FUNCTION:         DESCRIPTION:
#  .check.par		    Checks, whether the parameter are correctly specified, i.e., ascending ordered from lowest to the highest hierarchical level. (Internal function)
#  .check			      Checks, whether the parameters of two subsequent nodes are correctly specified. (Internal function)
#  hac.full			    Constructs 'hac' objects for fully nested Archimedean Copulae.
#  hac				      Constructs 'hac' objects for arbitrary nested Archimedean Copulae.
#  print.hac     	  Prints 'hac' objects.
#  hac2nacopula     Converts an 'hac' object into a 'nacopula' object.
#  .tree2nacList    Converts a 'tree' structure into a 'nacList'. (Internal function) 
#  nacopula2hac     Converts a 'nacopula' object into an 'hac' object.
#  .nacopula2tree   Converts a 'nacopula' object into a 'tree' structure. (Internal function)
##########################################################################################################################

.check.par = function(x){
	if((x$type == 2) || (x$type == 1) || (x$type == 8) || (x$type == 7)){
		ober.theta = 1}
	else
	if((x$type = 4) || (x$type = 3) || (x$type = 6) || (x$type = 5)){
		ober.theta = 1e-10}
	if((x$type = 10) || (x$type = 9)){
		ober.theta = 0}
	
	L = x$tree
	.check(L = L, theta = ober.theta)
}

#------------------------------------------------------------------------------------------------------------------------

.check = function(L, theta){
	n = length(L)
	if(L[[n]] < theta){
		return(warning("The dependency parameter of the nested AC should be higher than the parameter at the initial node."))
	}else{
	for(i in 1:(n-1)){
		if(class(L[[i]]) == "list"){
			.check(L = L[[i]], theta = L[[n]])
	}}}
}

#------------------------------------------------------------------------------------------------------------------------

hac.full = function(type, y, theta){
	n = length(y)
	if(n != (length(theta) + 1)){return(warning("The input arguments does not fit to a fully nested HAC"))}
	
	tree = list(y[n], y[n-1], theta[n-1])
	for(i in (n-2):1){
		tree = list(tree, y[i], theta[i])
	}
	hac(type = type, tree = tree)
}

#--------------------------------------------------------------------------------------------
 
hac = function(type, tree){
 	object = list(type = type, tree = tree)
 	class(object) = "hac"
 	if(((type == 10) || (type == 9)) && (max(.read.params(tree)) >= 1)){return(warning("The largest parameter of the Ali-Mikhail-Haq family must be < 1"))}
    .check.par(object);
 	object
}

#------------------------------------------------------------------------------------------------------------------------
	
print.hac = function(x, digits = 2, ...){
        cat("Class: hac", "\n", ...)
	if((x$type == 1) | (x$type == 2)){
 		   cat("Generator: Gumbel", "\n", ...)
 		   cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
	if((x$type == 3) | (x$type == 4)){
 		   cat("Generator: Clayton", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 	if((x$type == 5) | (x$type == 6)){
 		   cat("Generator: Frank", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 	if((x$type == 7) | (x$type == 8)){
 		   cat("Generator: Joe", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 	if((x$type == 9) | (x$type == 10)){
 		   cat("Generator: Ali-Mikhail-Haq", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 }

#------------------------------------------------------------------------------------------------------------------------

hac2nacopula = function(x){
	.family = character(1)
	.names = .get.leaves(x$tree)

	if((x$type == 1) | (x$type == 2)){
		.family = "G"
	}else
	if((x$type == 3) | (x$type == 4)){
		.family = "C"
	}else
	if((x$type == 5) | (x$type == 6)){
		.family = "F"
	}else
	if((x$type == 7) | (x$type == 8)){
		.family = "J"
	}else
	if((x$type == 9) | (x$type == 10)){
		.family = "A"
	}

    true.numbers = as.numeric(.names)
    if(any(is.na(true.numbers))){
       d = length(.names)
       for(j in 1:d){print(paste(.names[j], " <-> ", j, sep = ""))}
	   onacopulaL(.family, .tree2nacList(x$tree, .names, 1:length(.names)))
    }else{
       onacopulaL(.family, .tree2nacList(x$tree, .names, true.numbers))
    }
}

#------------------------------------------------------------------------------------------------------------------------

.tree2nacList = function(tree, names = NULL, numbers){	
	 if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
	 select = which(names %in% tree[-n])
	 
     if(any(s)){
         if(any(!s)){	
           res = list(tree[[n]], numbers[select], lapply(X = tree[which(!s)], FUN = .tree2nacList, names = names, numbers = numbers))
         }else{
           res = list(tree[[n]], numbers[select])
         }}else{
           res = list(tree[[n]], NULL, lapply(tree[which(!s)], FUN = .tree2nacList, names = names, numbers = numbers))
         }
     res
}

#------------------------------------------------------------------------------------------------------------------------

nacopula2hac = function(outer_nacopula){
  if(class(outer_nacopula) != "outer_nacopula"){stop("An outer_nacopula object is required.")}
	.family = character(1)
	n.childs = length(outer_nacopula@childCops)
	if(outer_nacopula@copula@name == "Gumbel"){
	   .family = if(n.childs > 0){1}else{2}
	}else
	if(outer_nacopula@copula@name == "Clayton"){
		.family = if(n.childs > 0){3}else{4}
  }else
	if(outer_nacopula@copula@name == "Frank"){
		.family = if(n.childs > 0){5}else{6}
	}else
	if(outer_nacopula@copula@name == "Joe"){
		.family = if(n.childs > 0){7}else{8}
	}else
	if(outer_nacopula@copula@name == "AMH"){
		.family = if(n.childs > 0){9}else{10}
	}
	hac(type = .family, tree = .nacopula2tree(outer_nacopula))
}

#------------------------------------------------------------------------------------------------------------------------

.nacopula2tree = function(outer_nacopula){
	 
   comps = outer_nacopula@comp
   n.comps = length(comps)
   n.childs = length(outer_nacopula@childCops)
	 
     if(n.childs > 0){
         if(n.comps > 0){	
           res = vector("list", n.comps + n.childs + 1)
           for(j in 1:n.comps){res[[j]] = paste(comps[j], sep = "")}
           for(j in (n.comps + 1):(n.comps + n.childs)){res[[j]] = .nacopula2tree(outer_nacopula@childCops[[j - n.comps]])}
           res[[n.comps + n.childs + 1]] = outer_nacopula@copula@theta
         }else{
           res = vector("list", n.childs + 1)
           for(j in 1:n.childs){res[[j]] = .nacopula2tree(outer_nacopula@childCops[[j]])}
           res[[n.childs + 1]] = outer_nacopula@copula@theta
         }}else{
           res = vector("list", n.comps + 1)
           for(j in 1:n.comps){res[[j]] = paste(comps[j], sep = "")}
           res[[n.comps + 1]] = outer_nacopula@copula@theta
         }
     res
}
