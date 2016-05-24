
Daim.control <- function(method="boot", number = 100, replace = TRUE, boot.size = 1,
k = 10, k.runs = 1, dependency = list(var = NULL, keep.id = FALSE))
{
	if(is.null(dependency[[1]]))
	dependency <- NULL
	if(method == "boot"){
		if(number < 0)
		stop("\n The value of 'number' supplied is < 0 ! \n")
		if(boot.size > 1)
		stop("\n The value of 'boot.size' supplied is > 1 ! \n")
		ans <- list(method=method, nboot=number, replace=replace,
					boot.size=boot.size, dependency=dependency)
	}
	else if(method == "cv"){
		replace <- FALSE
		if(k < 0)
		stop("\n The value of 'k' supplied is < 0 ! \n")
		if(k.runs < 1)
		stop("\n The value of 'k.runs' supplied is < 1 ! \n")
		ans <- list(method=method, k=k, k.runs=k.runs, dependency=dependency)
	}
	else if(method == "bcv"){
		if(number < 0)
		stop("\n The value of 'number' supplied is < 0 ! \n")
		if(boot.size > 1)
		stop("\n The value of 'boot.size' supplied is > 1 ! \n")
		if(k < 0)
		stop("\n The value of 'k' supplied is < 0 ! \n")
		ans <- list(method=method, nboot=number, replace=replace, k=k, 
					boot.size=boot.size, dependency=dependency)
	}	
	else{
		stop("\nThe value of 'method' must be one of 'boot', 'cv' or 'bcv'\n")
	}
	ans
}


#######################################
### Help function for Bootstrap k-fold CV on cluster
###

Daim.cluster.bcv <- function(n,formula,model,data,N,boot.size,replace,xval)
{
	ID.B <- sample(1:N,N*boot.size,replace=replace)
	ID.CV <- split(ID.B,ID.B)
	xgr <- 1:xval
	id <- sample(rep(xgr, length = length(ID.CV)), length(ID.CV))
	prob.oob.cv <- id.oob.cv <- lab.oob.cv <- NULL
	for(j in xgr){
		test.id <- id == j
		train <- data[unlist(ID.CV[!test.id],use.names = FALSE),]
		test <- data[unlist(ID.CV[test.id], use.names = FALSE),]
		prob.oob.cv <- c(prob.oob.cv,model(formula,train,test))
		id.oob.cv <- c(id.oob.cv,unlist(ID.CV[test.id], use.names = FALSE))
		lab.oob.cv <- c(lab.oob.cv,as.numeric(test[[1]])-1)
	}
	list(testind = id.oob.cv,
		 prob.oob = prob.oob.cv,
		 lab.oob = lab.oob.cv)
}


#######################################
### Help function for Bootstrap on cluster
###

Daim.cluster.boot <- function(n, formula, model, data, N, boot.size, replace)
{
	mylist <- sample(1:N, N*boot.size, replace=replace)
	train <- data[mylist,]
	testid <- unique(mylist)
	test <- data[-testid,]
	list(testind = (1:N)[-testid], inbaggind = mylist,
		 prob.oob = model(formula,train,test),
		 lab.oob = as.numeric(test[[1]])-1)
}



Daim.cluster.cv0 <- function(i,formula,model,data,N,xval)
{
	xgr <- 1:xval
	id <- sample(rep(xgr, length = N), N)
	prob.oob <- lab.oob <- testind <- vector(mode="list",length=xval)
	k <- 1
	for(j in xgr){
		test.id <- id == j
		train <- data[!test.id,]
		test <- data[test.id,]
		testind[[k]] <- which(test.id)
		prob.oob[[k]] <- model(formula,train,test)
		lab.oob[[k]] <- as.numeric(test[[1]]) - 1
		k <- k+1
	}
	list(testind = testind, prob.oob = prob.oob, lab.oob = lab.oob)
}




Daim.cluster.cv <- function(i,formula,model,data,N,xval,id)
{
	test.id <- id == i
	train <- data[!test.id,]
	test <- data[test.id,]
	list(testind = which(test.id),
		 prob.oob = model(formula,train,test),
		 lab.oob = as.numeric(test[[1]]) - 1)
}




Daim.boot.index <- function(N, boot.size, replace, dep.ind = NULL)
{
	ind <- sample(1L:N, N*boot.size, replace=replace)
	if(!is.null(dep.ind)){
		ind <- unlist(dep.ind[ind],use.names=FALSE)
	}
	ind
}



Daim.depend <- function(x, dep)
{
	id.var <- which(names(x) == dep[[1]])
	dep.ind <- split(1L:length(x[[id.var]]), x[[id.var]])
	if(length(dep) > 1 && !dep[[2]])
		x <- x[,-id.var]
	N <- length(dep.ind)
	list(x,N,dep.ind)
}




