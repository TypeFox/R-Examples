# functions.r ############################################################################################################
# FUNCTION:         DESCRIPTION:
#  tau2theta       	Transforms Kendall's rank correlation coefficient to the dependence parameter of an Archimedean copula.
#  theta2tau        Transfrorms the dependence parameter of an Archimedean copula to Kendall's rank correlation coefficient.
#  phi        		The generator function of Archimedean copula.
#  phi.inv			The inverse of the generator function.
#  copMult        	Computes the value of d-dimensional AC for a given sample with values in [0,1]^d.
#  par.pairs        Returns the pairwise arranged parameter in a matrix, so that the parameters correspond to the lowest hierarchical level at which the variables are joined. 
#  .pair.matr       Supplementary function of par.pairs. Arranges the variables pairwise and returns the corresponding value. (Internal function)
##########################################################################################################################

tau2theta = function(tau, type){
	n = length(tau)
	for(i in 1 : n){if((tau[i] < 0) | (tau[i] > 1)){return(warning(paste("tau[", i,"] should be in [0, 1).")))}}
    if((type == 2) || (type == 1))
        copGumbel@iTau(tau)
    else if((type == 4) || (type == 3))
        copClayton@iTau(tau)
	else if((type == 6) || (type == 5))
        copFrank@iTau(tau)
    else if((type == 8) || (type == 7))
        copJoe@iTau(tau)
    else if((type == 10) || (type == 9))
        copAMH@iTau(tau)
}

#-------------------------------------------------------------------------------------------------------------------------------

theta2tau = function(theta, type){
	n = length(theta)
    if((type == 2) || (type == 1)){
		for(i in 1 : n){if(theta[i] < 1){return(warning(paste("theta[", i,"] >= 1 is required.")))}}
        copGumbel@tau(theta)
    }else if((type == 4) || (type == 3)){
	    for(i in 1 : n){if(theta[i] <= 0){return(warning(paste("theta[", i,"] > 0 is required.")))}}
        copClayton@tau(theta)
    }else if((type == 6) || (type == 5)){
	    for(i in 1 : n){if(theta[i] <= 0){return(warning(paste("theta[", i,"] > 0 is required.")))}}
        copFrank@tau(theta)
    }else if((type == 8) || (type == 7)){
	    for(i in 1 : n){if(theta[i] < 1){return(warning(paste("theta[", i,"] >= 1 is required.")))}}
        copJoe@tau(theta)
    }else if((type == 10) || (type == 9)){
	    for(i in 1 : n){if((theta[i] >= 1) || (theta[i] < 0)){return(warning(paste("theta[", i,"] >= 0 and < 1 is required.")))}}
        copAMH@tau(theta)
   }
}

#-------------------------------------------------------------------------------------------------------------------------------

phi = function(x, theta, type){
	n = length(x)
	for(i in 1:n){if(x[i] < 0){return(warning(paste("x[", i,"] >= 0 is required.")))}}
	
    if((type == 2) || (type == 1)){
		if(theta >= 1){
			exp(-x^(1 / theta))
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}
	}else
	if((type == 4) || (type == 3)){
		if(theta > 0){
			(x + 1)^(-1 / theta)
		}else{
			return(warning(paste("theta > 0 is required.")))
		}
	}else
	if((type == 6) || (type == 5)){
		if(theta > 0){
			-log(-expm1(-x) + exp(-theta-x))/theta
		}else{
			return(warning(paste("theta > 0 is required.")))
		}			
	}else
	if((type == 8) || (type == 7)){
		if(theta >= 1){
			1 - (-expm1(-x))^(1 / theta)
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}			
	}else
	if((type == 10) || (type == 9)){
		if((theta >= 0) && (theta < 1)){
			(1 - theta) / (exp(x) - theta)
		}else{
			return(warning(paste("theta >= 0 and < 1 is required.")))
		}			
	}
}

#-------------------------------------------------------------------------------------------------------------------------------

phi.inv = function(x, theta, type){
	n = length(x)
	for(i in 1 : n){if((x[i] < 0) | (x[i] > 1)){return(warning(paste("x[", i,"] >= 0 and =< 1 is required.")))}}
	
    if((type == 2) || (type == 1)){
		if(theta >= 1){
			(-log(x))^theta
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}
	}else
	if((type == 4) || (type == 3)){
		if(theta > 0){
			(x^(-theta) - 1)
		}else{
			return(warning(paste("theta > 0 is required.")))
		}
	}else
	if((type == 6) || (type == 5)){
		if(theta > 0){
			 -log1p(exp(-theta)*expm1(theta-x*theta)/expm1(-theta))
		}else{
			return(warning(paste("theta > 0 is required.")))
		}
	}else
	if((type == 8) || (type == 7)){
		if(theta >= 1){
			-log1p(- (1 - x)^theta)
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}
	}else
	if((type == 10) || (type == 9)){
		if((theta >= 0) && (theta < 1)){
			log((1 - theta)/x + theta)
		}else{
			return(warning(paste("theta >= 0 and < 1 is required.")))
		}			
	}
}

#-------------------------------------------------------------------------------------------------------------------------------

copMult = function(X, theta, type){	
	phi(rowSums(phi.inv(X, theta = theta, type)), theta = theta, type)
}

#------------------------------------------------------------------------------------------------------------------------

par.pairs = function(hac, FUN = NULL, ...){
    tree = hac$tree
    vars = .get.leaves(tree)
    d = length(vars)
    matr = matrix(NA,nrow=d,ncol=d); colnames(matr)=rownames(matr)=vars
    matr = .pairs.matr(tree, matr)
    diag(matr) = NA
    
    if(class(FUN)=="function"){
        matr[lower.tri(matr)] = FUN(matr[lower.tri(matr)], ...)
        matr[upper.tri(matr)] = FUN(matr[upper.tri(matr)], ...)
    }else{
        if(!is.null(FUN)){
           if(FUN == "TAU")
            matr[lower.tri(matr)] = theta2tau(matr[lower.tri(matr)], type=hac$type)
            matr[upper.tri(matr)] = theta2tau(matr[upper.tri(matr)], type=hac$type)
        }}
    matr
}

#------------------------------------------------------------------------------------------------------------------------

.pairs.matr = function(tree, matr){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
     
     if(any(s)){
         if(any(!s)){
            l = sapply(tree[which(!s)], .get.leaves)
            for(i in 1:(length(l)-1))for(j in (i+1):length(l))matr[unlist(l[[i]]),unlist(l[[j]])]=matr[unlist(l[[j]]),unlist(l[[i]])]=tree[[n]]
            for(i in 1:length(l))matr[unlist(l[[i]]),unlist(tree[which(s)])]=matr[unlist(tree[which(s)]),unlist(l[[i]])]=tree[[n]]
            matr[unlist(tree[which(s)]), unlist(tree[which(s)])]=tree[[n]]
            for(i in which(!s)){
                matr = .pairs.matr(tree[i], matr)         
            }
         }else{
            matr[unlist(tree[-n]), unlist(tree[-n])]=tree[[n]]
         }}else{
            l = sapply(tree[-n], .get.leaves)
            for(i in 1:(length(l)-1))for(j in (i+1):length(l))matr[unlist(l[[i]]),unlist(l[[j]])]=matr[unlist(l[[j]]),unlist(l[[i]])]=tree[[n]]
            for(i in 1:(n-1)){
                 matr = .pairs.matr(tree[i], matr)            
            }
        }
    return(matr)
}