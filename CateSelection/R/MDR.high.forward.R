MDR.high.forward <-
function(x,y,order=NULL,trace=NULL,alpha=NULL,beta=NULL,pvalue=NULL,r2=NULL,...){

	if(is.null(trace))trace <- FALSE
	if(is.null(alpha)) alpha <- 0.1
	if(is.null(beta)) beta <- 0.05
	if(is.null(pvalue)) pvalue <- 0.01
	if(is.null(r2)) r2 <- 0.02
	if(is.null(order)) order <- 3
	res <- MDR.high.forward1(x=x, y=y, order=order, trace=trace, alpha=alpha, beta=beta, pvalue=pvalue, r2=r2)
	return(res)

}


MDR.high.forward1 <-
function(x,y,order,trace,alpha,beta,pvalue,r2,...){
	n <- ncol(x)
	index <- t(combn(n,order))
	if(order<3)stop("The order of interaction should be greater than 3")
	c <- order

	## Stage (1) MDR and MLR comparsion
	res1 <- MDR.sing.mod(x,y,order,trace)
	res1 <- res1[,-(1:order)]
	D1 <- res1[,1] - res1[,3]
	id <- which(D1 >= alpha)
	S1 <- matrix(index[id,],,c) ## new index
	MDR.R2 <- res1[id,1] ## MDR r-squared
	index <- NULL
	res1 <- NULL

	## Stage (2) Lower order interactions
	res2 <- low.anova(x,y,index=S1,trace)
	D2 <- MDR.R2 - res2[,1]
	id <- which(D2 >= beta)
	S2 <- matrix(S1[id,],,c)
	res2 <- NULL
	
	## Stage (3) Forward Selection
	res3 <- MDR_forward(Index=S2, dat=data.frame(x,y=y),alpha=pvalue,rsquared=r2)
	
	res <- unlist(res3[[1]])[,1:(c+1)]
	rownames(res) <- NULL
	colnames(res) <- NULL
	res <- matrix(res,,(c+1))	
	RES <- list(index = res[,1:c], R2 = res[,c+1])
	return(RES)

}
