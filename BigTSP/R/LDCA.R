library(glmnet)
library(tree)
library(randomForest)
library(gbm)
#input X, convert it to tsp
tspfeature = function(X)
{
	n=dim(X)[1]
	d=dim(X)[2]
	newX=NULL
	for(i in 1:(d-1))
	{
		for(j in (i+1):d)
		{
			newX=cbind(newX,as.numeric(X[,i]<X[,j]))
		}
	}
	#newX=cbind(1,newX)
	return(newX)
}

#LDCA <- function(x,...) UseMethod("LDCA")

LDCAEst <-function(X,y,nlambda1=100,lambda1=NULL,threshold=1e-07)
{
	tspf=tspfeature(X)
	fit=glmnet(tspf,y,family="binomial",alpha=1,nlambda=nlambda1,lambda=lambda1,maxit=10^5,thresh=threshold)
	i0=fit$a0
	beta=fit$beta
	lambda=fit$lambda
	df=fit$df
	return(list("a0"=fit$a0,"beta"=fit$beta,"lambda"=fit$lambda,"df"=fit$df,"dim"=fit$dim,"dev.ratio"=fit$dev.ratio,"nulldev"=fit$nulldev,"npasses"=fit$npasses,"jerr"=fit$jerr,"offset"=fit$offset,"classnames"=fit$classnames,"nobs"=fit$nobs))
}

#the main LDCA function
LDCA = function(X,y,nlambda=100,lambda=NULL,threshold=1e-07)
{
	Est = LDCAEst(X,y,nlambda,lambda,threshold)
	Est$call = match.call()
	class(Est) = c("LDCA","glmnet")
	Est
}

#cv.LDCA = function(x,...) UseMethod("cv.LDCA")

cv.LDCA = function(X,y,lambda=NULL,nfolds)
{
	tspf=tspfeature(X)
	cvfit=cv.glmnet(tspf,y,weights=rep(1,dim(tspf)[1]),lambda=lambda,nfolds=nfolds)
	fit=cvfit$glmnet.fit
	class(fit)="LDCA"
	cvobj=list("lambda"=cvfit$lambda,"cvm"=cvfit$cvm,"cvsd"=cvfit$cvsd,"cvup"=cvfit$cvup,"cvlo"=cvfit$cvlo,"nzero"=cvfit$nzero,"name"=cvfit$name,"LDCA.fit"=fit,"lambda.min"=cvfit$lambda.min,"lambda.lse"=cvfit$lambda.lse)
	cvobj$call = match.call()
	class(cvobj) = c("cv.LDCA","cv.glmnet")
	cvobj
}

print.LDCA = function(x,...)
{
	cat("Call:\n")
	print(x$call)
	cat("Lambda values:\n")
	print(x$lambda)
	cat("\nCoefficients:\n")
	print(x$beta)
}

print.cv.LDCA = function(x,...)
{
	cat("Call:\n")
	print(x$call)
	cat("The lambda value that gives the minimum cvm:\n")
	print(x$lambda.min)
	cat("The number of nonzero coefficients:\n")
	print(x$nzero)
}

predict.LDCA = function(object,newx,s=NULL,type=c("link","response","coefficients","nonzero","class"),exact=FALSE,offset,...)
{
	newx=tspfeature(newx)
	class(object)="glmnet"
	predict(object,newx,s=s,type=type,exact=exact,offset=offset,...)
}  

predict.cv.LDCA=function(object,newx,s=c("lambda.lse","lambda.min"),...){
  if(is.numeric(s))lambda=s
  else 
    if(is.character(s)){
      s=match.arg(s)
      lambda=object[[s]]
    }
    else stop("Invalid form for s")
    #newx=tspfeature(newx)
  predict(object$LDCA.fit,newx,s=lambda,...)
}

#tsp.tree = function(x,...) UseMethod("tsp.tree")
#start tree
tsp.tree = function(X,response,control=tree.control(dim(X)[1],...),method="recursive.partition",split=c("deviance", "gini"), x = FALSE, y = TRUE, wts = TRUE, ...)
{
	newx=tspfeature(X)
	d = data.frame(response,newx)
	rt=tree(response~.,data=d,control=control,method=method,split=split,x=x,y=y,wts=wts,...)
	rt$call=match.call()
	rt$data=d
	class(rt)=c("tsp.tree","tree")
	rt
}

predict.tsp.tree = function(object,newdata,type=c("vector", "tree", "class", "where"),split = FALSE, nwts, eps = 1e-3, ...)
{
	X=tspfeature(newdata[,2:dim(newdata)[2]])
	response=newdata[,1]
	newd=data.frame(response,X)
	class(object)="tree"
	predict(object, newd,type = type, split = split, nwts, eps = eps, ...)	
}

#start random forest

#tsp.randomForest = function(x,...) UseMethod("tsp.randomForest")
tsp.randomForest = function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree=500,
             type="classification",mtry=if (!is.null(y) && !is.factor(y))
             max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             replace=TRUE, classwt=NULL, cutoff, strata,
             sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
             maxnodes=NULL,
             importance=FALSE, localImp=FALSE, nPerm=1,
             proximity=FALSE, oob.prox=proximity,
             norm.votes=TRUE, do.trace=FALSE,
             keep.forest=!is.null(y) && is.null(xtest),
             keep.inbag=FALSE, ...) 
{
	newx=tspfeature(x)
	if(type=="classification") {y=as.factor(y)}
	else {y=y}
	rf=randomForest(newx, y=y,  xtest=xtest, ytest=ytest, ntree=ntree,
             mtry=mtry,
             replace=replace, classwt=classwt, cutoff, strata,
             sampsize = sampsize,
             nodesize = nodesize,
             maxnodes=maxnodes,
             importance=importance, localImp=localImp, nPerm=nPerm,
             proximity=proximity, oob.prox=oob.prox,
             norm.votes=norm.votes, do.trace=do.trace,
             keep.forest=keep.forest,
             keep.inbag=keep.inbag, ...)
    rf$call=match.call()
    class(rf)=c("tsp.randomForest","randomForest")
    rf
}

predict.tsp.randomForest = function(object,newdata,type="response",norm.votes=TRUE,predict.all=FALSE,proximity=FALSE,nodes=FALSE,cutoff,...)
{
	X=tspfeature(newdata)
	class(object)="randomForest"
	predict(object, X, type="response",
  norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE,
  cutoff, ...)
}

#start gradient boosting method

#tsp.gbm = function(x,...) UseMethod("tsp.gbm")

tsp.gbm = function(x,y,offset=NULL,misc=NULL,distribution="bernoulli",w=NULL,var.monotone=NULL,n.trees=100,interaction.depth=1,n.minobsinnode=10,shrinkage=0.001,bag.fraction=0.5,train.fraction=1.0,keep.data=TRUE,verbose=TRUE)
{
	newx=tspfeature(x)
	rg=gbm.fit(newx,y,offset=offset,misc=misc,distribution=distribution,w=w,var.monotone=var.monotone,n.trees=n.trees,interaction.depth=interaction.depth,n.minobsinnode=n.minobsinnode,shrinkage=shrinkage,bag.fraction=bag.fraction,train.fraction=train.fraction,keep.data=keep.data,verbose=verbose)	
	rg$call=match.call()
	class(rg)=c("tsp.gbm","gbm")
	rg
}
predict.tsp.gbm = function(object,newdata,n.trees,type="link",single.tree=FALSE,...)
{
	newx=tspfeature(newdata)
	class(object)="gbm"
	predict(object,newx,n.trees=n.trees,type=type,single.tree=single.tree,...)
}
