"trudi" <-
function(mat, q = 0, weights = TRUE, base = exp(1)){
	## gamma pi
	pi.g <- colSums(mat/sum(mat))
	## alternative gamma pi
	#pi.ga <- mat/sum(mat)
	## alpha pi
	pi.a <- mat/rowSums(mat)
	## weights (with weights = TRUE and q = 0 its still diversity
	## as opposed to with q = 0 and weights = FALSE its richness based)
	if(weights){
		su.w <- rowSums(mat/sum(mat))
		if(q==1){
			gm.pi.a <- apply(pi.a, 1, function(x) sum(x*log(x, base=base), na.rm=TRUE))
			td.a <- exp(sum(-su.w * gm.pi.a))
		}
		else{
			## true alpha
			#gm.pi.aa <- apply(pi.a, 1, function(x) (sum(x*(x^(q-1)), na.rm=TRUE))^(1/(q-1)))
			gm.pi.a <- apply(pi.a, 1, function(x) (sum(x*(x^(q-1)), na.rm=TRUE)))
			gm.pi.a <- sum(gm.pi.a * su.w)^(1/(q-1))
			td.a <- 1/gm.pi.a
		}
	}
	else{
		if(q!=0){warning("giving no weights in case q != 0")}
		## true alpha (richness case)
		td.a <- mean(rowSums(mat))
		if(q==1){
			gm.pi.a <- apply(pi.a, 1, function(x) sum(x*log(x, base=base), na.rm=TRUE))
			td.a <- exp(-sum(gm.pi.a/nrow(mat)))
		}
	}
	## true gamma
	if(q==1){
		#td.g <- sum(apply(pi.ga, 2, function(x) -sum(x*su.w) * log(sum(x*su.w), base=base)))
		td.g <- exp(-sum(pi.g*log(pi.g, base=base), na.rm=TRUE))
	}
	else{
		gm.pi.g <- sum(pi.g*(pi.g^(q-1)), na.rm=TRUE)^(1/(q-1))
		td.g <- 1/gm.pi.g
	}	
	## true beta
	td.b <- td.g/td.a
	return(unlist(list(gamma=td.g, beta=td.b, alpha=td.a)))
}