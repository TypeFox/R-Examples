ols.refit <-
function(y, X, mainind, interNames, family){
	Xall = X[,mainind]
    if(length(interNames)>0){
    #	browser()
	for(i in 1:length(interNames)){
	pair = as.numeric(strsplit(interNames[i], "X")[[1]][2:3])
	xi = X[,pair[1]]*X[,pair[2]]
	Xall = cbind(Xall,xi)		
	}
	}
	#browser()
	lmfit = glm(y~Xall,family=family)
	beta = coef(lmfit)
	link=as.vector(cbind(1,Xall)%*%beta)
	n = length(y)
	if(family=='poisson') loglik = -2*sum(exp(link)+2*y*link)
	if(family=='binomial') loglik = 2*sum(log(1+exp(link))-y*link)
    if(family=='gaussian') loglik = n*log(mean((y-link)^2))
    #loglik = deviance(lmfit)
	return(list(beta=coef(lmfit), loglik=loglik, link=link))
}
