sing.mod <-
function(x,y,order=NULL,alpha=NULL,beta=NULL,delete=NULL,trace=NULL,...){
	if(is.null(alpha))alpha <- 0.01 # p-value
	if(is.null(beta))beta <- 0.01 # R-squared
	if(is.null(delete))delete <- FALSE
	if(is.null(trace))trace <- FALSE
	if(is.null(order)) order <- 2
	res <- sing.mod1(x=x,y=y,order=order,alpha=alpha,beta=beta,delete=delete,trace=trace)
	return(res)
}



sing.mod1 <-
function(x,y,order,alpha,beta,delete,trace){
	
	n <- ncol(x)
	index <- t(combn(n,order))

	##change column names
	n <- dim(x)[2]
	ColNames <- NULL
  	for(i in 1:n){ColNames[i] <- paste("x",i,sep="")}
	colnames(x) <- ColNames

	c <- ncol(index)
	r <- nrow(index)   
	y1 <- as.matrix(y,,1)
	dat <- x
	Pvalue <- numeric()
	R <- numeric()
	R_Adjust <- numeric()

	for(i in 1:r){
		if(trace)cat("Step:",i,"/",r,"\n")
		if(c<1)stop("'c' must be c >= 1")
    		if(c==1){
			ok <- complete.cases(dat[,index[i]])
			X <- factor(dat[,index[i]])
    		}
		if(c>1){
			ok <- complete.cases(dat[,index[i,]])
			X <- dat[,index[i,1]]
			for(j in 1:(c-1))X <- paste(X,dat[,index[i,j+1]],sep=":")	    
			X <- factor(X)	  
		}

		Y <- y1[ok,]
		X <- X[ok]

		reg <- lm(Y~X)
		B <- summary(reg)
		r2 <- B$r.squared
		R[i] <- r2
		nR <- length(Y)
		pR <- 1
		R_Adjust[i] <- 1-(1-r2)*(nR-1)/(nR-pR-1)
		F <- B[[10]][[1]]
		df1 <- B[[10]][[2]]
		df2 <- B[[10]][[3]] 
		Pvalue[i] <- 1-pf(F,df1,df2)
	}

	NewInd <- cbind(index,Pvalue,R,R_Adjust)
	if(delete){
		id1 <- which(Pvalue>alpha)
		id2 <- which(R<beta)
		id0 <- union(id1,id2)
		if(length(id0)>0)NewInd <- NewInd[-id0,]
	}
	return(NewInd)
}
