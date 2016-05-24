getResponseCat1 <- function(mat.prob, n.cat, n.obs, beta0, n.ias, sample.y=TRUE){
	mat.prob <- exp(mat.prob) / (1 + exp(mat.prob))
	mat.ias <- matrix(0, n.obs, n.cat)
	for(i in 1:n.cat)
		mat.ias[,i] <- mat.prob[,i] > exp(beta0[i]) / (1 + exp(beta0[i]))
	cl <- vec.which <- numeric(n.obs)
	ids.none <- rowSums(mat.ias) == 0
	cl[ids.none] <- sample(n.cat,sum(ids.none), replace=TRUE, prob=exp(beta0)/(1+exp(beta0)))
	if(!sample.y)
		cl[!ids.none] <- vec.which[!ids.none] <- max.col(mat.prob[!ids.none,],
			ties.method="random")
	else{
		for(i in 1:n.cat){
			ids <- which(mat.ias[,i]==1)	
			vec.which[ids] <- i
			tmp <- mat.prob[ids,]
			tmp[,-i] <- tmp[,-i] / rowSums(tmp[,-i]) * (1 - tmp[,i])
			if(n.ias[i] == 1)
				cl[ids] <- sample(n.cat, length(ids), replace=TRUE, prob=tmp[1,])
			else{
				uni.prob <- unique(tmp[,i])
				tmp.cl <- numeric(length(ids))
				for(j in 1:length(uni.prob)){
					ids2 <- which(tmp[,i] == uni.prob[j])
					tmp.cl[ids2] <- sample(n.cat, length(ids2), replace=TRUE,
						prob=tmp[ids2[1],])
				}
				cl[ids] <- tmp.cl
			}
		}
	}
	return(list(cl=cl, vec.which=vec.which))
}


getResponseRef <- function(mat.prob, n.cat, n.obs, beta0, sample.y=TRUE){
	mat.ias <- matrix(0, n.obs, n.cat)
	for(i in 1:n.cat)
		mat.ias[,i] <- mat.prob[,i] > beta0[i]
	ids <- rowSums(mat.ias) > 0
	vec.which <- numeric(n.obs)
	vec.which[ids] <- max.col(mat.prob[ids,])
	mat.prob2 <- cbind(1, exp(mat.prob))
	mat.prob2 <- mat.prob2 / rowSums(mat.prob2)
	if(!sample.y)
		cl <- max.col(mat.prob2, ties.method="random") - 1
	else
		cl <- apply(mat.prob2, 1, function(x) sample(0:n.cat, 1, prob=x))
	return(list(cl=cl, vec.which=vec.which))
}






