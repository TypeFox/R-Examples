scair <- function(x,y,shape=rep("l",1), family=gaussian(), weights=rep(1,length(y)), epsilon=1e-8, delta=0.1, indexgen=c("unif", "norm"), iter = 200, allnonneg = FALSE){

     ## Some checking on inputs
        if(is.vector(x)==TRUE) x=matrix(x,ncol=1)
        y = as.vector(y)
        shape = as.vector(shape)
	if(is.matrix(x)==FALSE||is.numeric(x)==FALSE||is.vector(y)==FALSE||is.numeric(y)==FALSE) stop("Input error!")
        if((length(y)!=nrow(x))||(length(weights)!=nrow(x))) stop("Input error: x, y, d mismatch.")
        d=NCOL(x)
	m=length(shape)
        if(m > d) stop("Input error: x and shape mismatch.")
        N=NROW(x)
        if(sum(shape=="l") > 1) stop("At most one linear component is allowed for identifiability.")
	if(delta >= 1) stop("Reasonable delta, if required, should usually be between 0.01 and 0.95.")
	indexgen <- match.arg(indexgen)

     ## Generate index matrices and compute log-likelihood
	mindeviance = Inf
	index = matrix(0,ncol=m,nrow=d)
	count = 0
	while(count < iter){
		if(indexgen == "unif") w = matrix(runif(d*m)*2-1,ncol=m) else w = matrix(rnorm(d*m),ncol=m)
		if(allnonneg == TRUE) w = abs(w)
		# normalisation
		linearcomponent = 0
		for (j in 1:m) if(shape[j]=="l") linearcomponent = j
		if (linearcomponent!=0) {for(j in (1:m)[-linearcomponent]) w[,j] = w[,j] - w[,linearcomponent]*sum(w[,j]*w[,linearcomponent])/sum(w[,linearcomponent]^2)}
		for (j in 1:m) w[,j] = w[,j]/sum(abs(w[,j]))
		for (j in 1:m) if(shape[j]=="l"||shape[j]=="cvx"||shape[j]=="ccv") {if (w[1,j]<0) w[,j] = -w[,j]}

     ## try to determine whether delta is needed here - check whether the shape constraints imply a convex / concave function
		ccvccx = 0		
		if (prod((shape == "cvx") + (shape == "cvxin") + (shape == "cvxde") + (shape == "l")) == 1 ) ccvccx = 1
		if (prod((shape == "ccv") + (shape == "ccvin") + (shape == "ccvde") + (shape == "l")) == 1 ) ccvccx = 1
     ## and check whether the ridge functions are either all increasing or all decreasing
                inde = 0
		if (prod((shape == "in") + (shape == "cvxin") + (shape == "ccvin")) == 1 ) inde = 1
		if (prod((shape == "de") + (shape == "cvxde") + (shape == "ccvde")) == 1 ) inde = 1
		if((m == 1) || (ccvccx == 1) || (inde == 1 && allnonneg == TRUE) || (min(eigen(t(w)%*%(w))$values) >= delta) ){
			output = scar(x %*% w,y,shape=shape,family=family,weights=weights,epsilon=epsilon)
			if(count == 0 || (output$deviance < mindeviance)) {
				obj = output
				index = w
				mindeviance = output$deviance
			} 
			count = count + 1			
		} 
	}
	result = list(x=x,y=y,index = index, shape=shape, weights=weights, family=family, componentfit=obj$componentfit, constant=obj$constant, deviance=obj$deviance, nulldeviance=obj$nulldeviance, delta = delta, iter=iter, allnonneg=allnonneg)
	class(result) <- "scair"
        return(result)
}
