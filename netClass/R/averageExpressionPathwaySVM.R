#' Cross validation for aepSVM (aepSVM) using SAM to select significant differential expressed genes
#'
#' Implementation of the Aeverage Expression of Pathway( aepSVM) algorithm.

#'
#' @param x: a p x n matrix of expression measurements with p samples and n genes.
#' @param y: a factor of length p comprising the class labels.
#' @param DEBUG: show debugging information in screen or not.
#' @param scale: a character vector defining if the data should be centered and/or scaled.
#' Possible values: are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
#' @param Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @param Gsub: an adjacency matrix that represents the underlying biological network.
#' @seed  seed: for random sampling.
#' @return a list with the results of the cross-validation.  \item{feat}{the selected features}  \item{fit}{the fitted SVM model} \item{auc}{the prediction auc of the fitted SVM model} \item{lables}{the test labels}  
#' @references Guo Z, et al.: Towards precise classification of cancers based on robust gene functional expression profiles. \emph{BMC Bioinformatics} 2005, 6:58.
#' @export
#' @note The optimal number of features is found by using the span estimate. See Chapelle, O., Vapnik, V., Bousquet, O., and Mukherjee, S. (2002). Choosing multiple parameters for support vector machines. \emph{Machine Learning}, 46(1), 131-159.
#'
#' @author  Yupeng Cun \email{yupeng.cun@gmail.com}
#' @examples
#' library(netClass)
#' data(expr)
#' data(ad.matrix)
#' x <- expr$genes
#' y <- expr$y
#'
#' library(KEGG.db)
#' r.aep <- cv.aep(x, y, folds=5, repeats=3, parallel=FALSE,cores=2,Gsub=ad.matrix,Cs=10^(-3:3),seed=1234,DEBUG=TRUE)

cv.aep <- function(x, y, folds=10, repeats=5, parallel = FALSE, cores = 2, DEBUG=TRUE, Gsub=matrix(1,100,100), Cs=10^(-3:3), seed=1234)
{
	multicore <- ("package:parallel" %in% search())
  
	if(multicore == TRUE && parallel == TRUE)  
	{
    
		if(is.null(cores))       
			cores <- 2  	  		  
		parallel <- TRUE  
	}  
	else  
	{
		if(parallel == TRUE)       #
		cat('You nedd Package \'parallel\'\n')
		cat('Computation will performed sequential crossvalidation.\n', sep='')
		parallel <- FALSE
	}

	if(!is.factor(y)) stop("y must be factor!\n")
  
	if(length(levels(y)) != 2) stop('y must be factor with 2 classes.\n')
  
	if(length(y) != nrow(x)) stop('y must have same length to nrow(x).\n')
	
	
	int= intersect(colnames(x),colnames(Gsub))
  
	n     <- length(y)
	folds <- trunc(folds)

	if (folds < 2) stop("folds should be >= 2.\n")
	if (folds > n) stop("folds should be =< the number of sample size.\n")

	cuts  <- cv.repeats <- list()  	
	set.seed(seed)
	for(r in 1:repeats)
	{t
		perm <- sample(1:n) #Sampling a random integer between 1:n
		repeat.models <- NULL 		
    
		for(k in 1:folds) #randomly divide the training set in to 10 folds
		{
			tst <- perm[seq(k, n, by=folds)]  #      
			trn <- setdiff(1:n, tst)            
			cuts[[k]] <- list(trn=trn, tst=tst)    
		}    
		
		#pb <- txtProgressBar(min = 0, max = folds, style = 3)   #txtProgressBar id funtion in packages "utils", Text progress bar in the R console
    
		if(DEBUG) cat('Starting classification of repeat:',r,'\n')
    
		if(parallel)	repeat.models <- mclapply(1:folds, classify.aep , cuts=cuts, x=x, y=y, cv.repeat=r, int=int, DEBUG=DEBUG,Cs=Cs)
		else		repeat.models <-   lapply(1:folds, classify.aep , cuts=cuts, x=x, y=y, cv.repeat=r, int=int, DEBUG=DEBUG,Cs=Cs)
	
		#close(pb)
    
		if(length(repeat.models) != folds)
		{      
			geterrmessage()      
			stop("One or more processes did not return. May be due to lack of memory.\n")
		}
	
		if(DEBUG)	cat('All models of repeat:',r,'have been trained.\n')
    
		cv.repeats[[r]] <- repeat.models
  
  	} 
  	
  	
  	auc <- sapply(cv.repeats, function(cv.repeat) sapply(cv.repeat, function(model) model$auc))  
  	colnames(auc) <- paste("Repeat",1:repeats,sep="")  
  	rownames(auc) <- paste("Fold",1:folds,sep="")  
  	
  	fits <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$model))  
  	names(fits) <- paste("Repeat",1:repeats,sep="")  
  	fits <- lapply(fits, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })  
  	
  	genes <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$genes))  
  	names(genes) <- paste("Repeat",1:repeats,sep="")  
  	genes <- lapply(genes, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x }) 
  	
  	res <- list(auc=auc, fits=fits, feat= genes, labels=y)  
  	class(res) <- 'pathClassRes' 
  	return(res)  
	#return ( cv.repeats)
}

# Training and predicting using aepSVM (aepSVM) classification methods
#-----------------------------  
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels.
# DEBUG: show debugging information in screen or not.
# scale: a character vector defining if the data should be centered and/or scaled. 
#		Possible values: are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
# Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
# Gsub: an adjacency matrix that represents the underlying biological network.
#-----------------------------  
# Return a list with the results of the training model and predcit AUC   

classify.aep <- function(fold, cuts, Cs, x=x, y=y, cv.repeat, int=int, DEBUG=DEBUG,  Gsub=Gsub)
{
	gc()  	
	if(DEBUG) cat('starting Fold:',fold,'\n')  
	
	trn <- cuts[[fold]]$trn
	tst <- cuts[[fold]]$tst	

	train	<- train.aep(x=x[trn,], y=y[trn], DEBUG=DEBUG, int=int, Gsub=Gsub, Cs=Cs)
	test	<- predictAep(train= train, x=x[tst,], y=y[tst], DEBUG=DEBUG,  Gsub=Gsub)
		
	auc	= test$auc
	cv	= test$cv
	genes = train$sig.genes
	if(DEBUG)
	{ cat(" Test AUC =", auc, "\n")
		cat('Finished fold:',fold,'\n\n')
	}
	
	gc()
	list(fold=fold, model=train, auc=auc, cv=cv,genes=train$sig.genes)
}

train.aep <- function(x=x, y=y, DEBUG=FALSE, int=int, Gsub=Gsub, Cs=10^(-3:3))
{	
	#if(is.null(int)) int= intersect(colnames(x), colnames(Gsub))
	# select significant differential expressed genes
	#x=x[,int]
	#Gsub=Gsub[int,int]
	sigGens= selecteGenes(x=x, y=y, DEBUG=DEBUG)
	
	if(DEBUG) cat("There are " ,length(sigGens), " DE genes selected\n")

	# mean expression gene matrix of pathways
	trainPath <- probeset2pathway(x=x[,sigGens], int=int, sigGens=sigGens)	
	#print(ncol(trainPath))
	if(is.null(ncol(trainPath)))
	{
		sigGens   <- int
		trainPath  <- probeset2pathway(x=x[,sigGens], int=int,sigGens=sigGens)
	}
	
	# train the matrix of pathways in L2-SVM
	trained <- fit.aep( x=trainPath, y=y, DEBUG=DEBUG, scale="center", Cs=Cs)  
	# calculate the AUC	
	
	#trained$graph = graph.adjacency(Gsub[sigGens,sigGens], mode="undirected")	
	 
	res=list(trained=trained,sig.genes=sigGens,int=int)
	return(res)
}

predictAep <- function(train= train, x, y, DEBUG=FALSE,  Gsub=Gsub)
{
	label	<- sign(as.numeric(y) - 1.5) # because y is a factor  
	cv      <- double(length=length(y))
	
	sigGens	= train$sig.genes
	int		= train$int
	#x=x[,int]
	#Gsub=Gsub[int,int]
	# mean expression gene matrix of pathway
	testPath  <- probeset2pathway(x=x[,sigGens], int=int,sigGens=sigGens)	
	#predict the test matrix of pathways
	test    <- svm.predict(fit=train$trained$fit, newdata=testPath,type="response")#response		
	#cat("\n lable: ", label, "\n test: ", test,"\n")
	auc   <- calc.auc(test, label)	
	
	if(DEBUG) cat(" Test AUC =", auc, "\n")

	
	res=list(auc=auc, pred=test)
	return (res)
}




# Training trained fold using aepSVM methods
#-----------------------------  
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels.
# DEBUG: show debugging information in screen or not.
# scale: a character vector defining if the data should be centered and/or scaled. 
#		Possible values: are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
# Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
# Gsub: an adjacency matrix that represents the underlying biological network.
#-----------------------------  
# Return a list with the results of the training model for predciting   


fit.aep <- function(x, y, DEBUG=FALSE, scale=c('center', 'scale'),Cs=10^c(-3:3))
{ 
	best.bound = Inf
  
	feat = colnames(x)
	
	fit = svm.fit(x=x, y=y, Cs=Cs, scale="center", DEBUG=DEBUG)  
	          
	cat('Best Model is: Spanbound=',fit$error.bound,', C=',fit$C,',', length(fit$features),'features.\n')
  
	result <- list(features=feat, error.bound = fit$best.bound, fit=fit)
  	#class(result) <- 'sg'
	
	return(result)
}


### selected significant gene using SAM
selecteGenes <- function(x=x, y=y, DEBUG=DEBUG)
{ 
	best.bound = Inf	
  
	#change y labels form factor {-1,1} to numeric {1, 2}  
	yy <- sign(as.numeric(y)-1.5)    
	for(i in 1:length(yy))
		if(yy[i]==-1)yy[i]=0  
	yy=yy+1
	
	tx=t(x)
	# select significant genes here  
	d=list(x=tx,y=yy,genenames=as.character(rownames(tx)),logged2=TRUE)
	samr <- samr(d, resp.type="Two class unpaired", nperms=200)

	delta.table<-samr.compute.delta.table(samr)
	del<- delta.table[delta.table[,5] < 0.05, 1][1]
	if(is.na(del))
       del = max(delta.table[delta.table[,4] > 0, 1])
	#cat("\n del:",del,"\n")
	#cat("\nsigTable:\n==========================================\n",delta.table,"\n")
	
	siggenes.table<- samr.compute.siggenes.table(samr, del, d, delta.table)	
			
	genes.up = siggenes.table$genes.up[, 2]
	genes.down = siggenes.table$genes.lo[, 2]
   	   
	genes.sig = unique(c(genes.up, genes.down))  	
	cat("\n 1 st genes.sig:\t",length(genes.sig), "\t", "del:",del,"\n\n")
	
	if(length(genes.sig)<10)
	{
		del<- delta.table[delta.table[,5] < 0.5, 1][1]
		if(is.na(del))
        del = max(delta.table[delta.table[,4] > 0, 1])
	#cat("\n del:",del,"\n")
	#cat("\nsigTable:\n==========================================\n",delta.table,"\n")
		siggenes.table<- samr.compute.siggenes.table(samr, del, d, delta.table)	
			
		genes.up = siggenes.table$genes.up[, 2]
		genes.down = siggenes.table$genes.lo[, 2]
   	   
		genes.sig = unique(c(genes.up, genes.down))  	
		cat("\n 2ed genes.sig:\t",length(genes.sig), "\t", "del:",del,"\n\n")
		
	}
	
	pvSGids=genes.sig
	cat("\n pvSGids:\t",length(pvSGids), "\t", "del:",del,"\n\n")

	#if(length(pvSGids)!=0)pvSGids= siggenes.table$genes.up[1,2]
		
	feat = pvSGids	
	
	return(feat)
}

######
# generae a mean gene expression of genes of each pathway matrix
# ------------
# x: gene expression data
# int: common genes between pathway genes and genes in gene expression profile
# sigGens: significant gene expression using SAM methods
# ------------ 
probeset2pathway <-function(x=x, int=int,sigGens=sigGens)
{
	#
	# for (reactome.db)	
	#KEGGids = mget(int,reactomeKEGGEXTID2PATHID,ifnotfound=NA)	

	# for require(KEGG.db) ,AnnotationDbi::
	KEGGids =  mget(int,KEGGEXTID2PATHID,ifnotfound=NA)
	KEGGids = KEGGids[which(KEGGids!="NA")]

	##	
	kidList <- unlist(KEGGids)
	kidList <- kidList[which(kidList!="NA")]		
	ukidList <-unique(kidList)
	
	## form probset belond to KeegID, 	
	nCol = length(ukidList)
	nRow = nrow(x)
	nColx = ncol(x)	
	kse <- matrix(0, ncol=nCol, nrow=nRow)		
	i=1
	for(k in ukidList)	
	{		
		index = intersect(names(kidList[which(kidList==k)]), sigGens)				
		#cat("\nindex", index,"\n")		
		for(j in 1:nRow)
		{
			kse[j,i] = mean(x[j, index])
		}
		i=i+1
	}
	colnames(kse) = ukidList
	rownames(kse) <- rownames(x)
	cat("dim(kse):,\n",dim(kse),"\n")
	
	index=which(kse[2,]!="NaN")
	kse=kse[,index]
	cat("dim(kse):,\n",dim(kse),"\n")
	return (kse)
}



svm.fit = function(x, y, Cs, scale, DEBUG=FALSE){

  if(missing(x))     stop('No epxression matrix provided.')
  if(missing(y))     stop('No class-lables provided.')
  if(missing(Cs))    stop('No tuning parameter \'C\' provided.')
  if(missing(scale)) stop('Parameter \'scale\' must be in \'scale\', \'center\' or NULL.')
  if(length(levels(factor(y))) != 2) warning('y must have 2 levels.')
    
  scale.mean <- scale.std <- NULL

  if(!is.null(scale)){
    scale <- tolower(scale)

    if("center" %in% scale){
      x = scale(x,center=T)
      ## save centering coefficient
      scale.mean = attr(x,"scaled:center")
      names(scale.mean) = colnames(x)
    }
    
    if("scale" %in% scale){
      x = scale(x,center=F,scale=T)
      ## save scaling coefficient    
      scale.std = attr(x,"scaled:scale")
      names(scale.std) = colnames(x)
    }
  }

  ## this happens sometimes whenn all
  ## probe sets have exactely the same value
  ## and after centering everything is zero then.
  ## Due to this the scaling fails and produces
  ## NaNs:
  ## setting NaN columns to 0
  nan.cols     <- apply(x,2,function(y) all(is.nan(y)))
  x[,nan.cols] <- 0
  
  K <- kernelMatrix(vanilladot(), x)
  best.bound <- Inf

  for(C in Cs){
    if(DEBUG) cat('Trying C=',C,'\n')
    K2 = as.kernelMatrix(K + 1/C*diag(NROW(x)))
    fit.tmp = ksvm(K2, y, C=Inf, type="C-svc", shrinking=FALSE, tol=0.01,scaled=F)
    bound = spanbound(fit.tmp, K2, sign(as.numeric(y) - 1.5))
    if(bound < best.bound){
      model = fit.tmp
      best.bound = bound
      Cbest = C
    }
  }
  if(DEBUG) cat('Best C=',Cbest,'\n')
  
  ## alphaindex: The index of the resulting support vectors in the data matrix
  svs <- unlist(alphaindex(model))

  ## coef: The corresponding coefficients times the training labels.
  w <- abs(t(unlist(coef(model))) %*% x[svs,])
  
  fit <- list(fit.svm=model, w=w, K=K, C=Cbest, xsvs=x[svs,,drop=FALSE], error.bound=best.bound, scale.mean=scale.mean, scale.std=scale.std, features=colnames(x), R=NULL)	
  return(fit)
}

svm.predict = function(fit, newdata, type="response"){

  ## do the prediction only with those genes
  ## that were use for training
  newdata <- newdata[,fit$features]

  if(!is.null(fit$scale.mean))
    newdata <- scale(newdata, center=fit$scale.mean[fit$features], scale=FALSE)
  
  if(!is.null(fit$scale.std))
    newdata <- scale(newdata, center=FALSE, scale=fit$scale.std[fit$features])
  
  Ktst        <- kernelMatrix(vanilladot(), newdata, fit$xsvs[,fit$features])
  Ktst2       <- kernelMatrix(rbfdot(sigma=0.001), newdata, fit$xsvs[,fit$features])
  ident       <- which(Ktst2 == 1)		
  Ktst[ident] <- Ktst[ident] + 1/fit$C		
  alpha       <- as.matrix(unlist(coef(fit$fit.svm)))	
  yhat        <- Ktst%*%alpha - b(fit$fit.svm)

  if(type == "class") yhat <- sign(yhat)	

  return(yhat)
}


spanbound <- function(fit, xtrn, ytrn){
  
  svindex = unlist(alphaindex(fit))	
  alpha = unlist(coef(fit))	
  pos = which(alpha > 0)
  neg = which(alpha < 0)
  alpha = abs(alpha)	
  if(class(xtrn) != "kernelMatrix"){
    yhat = predict(fit, xtrn, type="decision")
    if(param(fit)$C != Inf)
      error("span bound is only for L2-SVM!")
    K = kernelMatrix(kernelf(fit), xtrn)
  }
  else{
    K = xtrn
    yhat = K[,svindex]%*%as.matrix(unlist(coef(fit))) - b(fit)
  }
  output = ytrn*yhat	
  Cpos = Inf
  Cneg = Inf
  eps = 1e-5	
  boundpos = (alpha[pos] >= Cpos*(1-eps))
  boundneg = (alpha[neg] >= Cneg*(1-eps))
  sv1pos = svindex[pos[!boundpos]]
  sv2pos = svindex[pos[boundpos]]
  sv1neg = svindex[neg[!boundneg]]
  sv2neg = svindex[neg[boundneg]]	
  sv1 = sort(c(sv1pos, sv1neg))
  sv2 = sort(c(sv2pos, sv2neg))
  n = ncol(K)
  span = double(n)
  alpha1 = double(n)
  alpha1[svindex] = alpha 		
  if(length(sv1) > 0){ # in-bound SVs 								
    ell = length(sv1)	
    invK = chol2inv(chol(K[sv1,sv1,drop=FALSE]))
    T = -1/sum(invK)
    T2 = invK%*%as.matrix(rep(1,ell))*T
    T3 = t(as.matrix(rep(1,ell)))%*%invK
    invKSV = rbind(cbind(invK + T2%*%T3, -T2), cbind(-T*T3, T))
    tmp = diag(as.matrix(invKSV)) + 1e-10
    span[sv1] = 1./tmp[1:ell]
  }	
  else
    warning("No in-bound SVs!")			
  if(length(sv2) > 0){	# bound SVs	
    span[sv2] = diag(as.matrix(K[sv2,sv2,drop=FALSE]))
    if(length(sv1) > 0){
      V = rbind(K[sv1,sv2,drop=FALSE], rep(1,length(sv2)))			
      span[sv2] = span[sv2] - diag(t(V)%*%invKSV%*%V)			
    }				 
  }	
  loo = mean((output - alpha1*span <= 0)*1)				
  ## cat("Span bound =", loo,"\n\n")		
  loo
}



calc.auc <- function(prob,labels)
{
  ## this corrects a bug in ROCR:
  ## if all labels are from one group and there
  ## is no missclassification, ROCR is not able
  ## to calculate the auc
  ## patch => add a artificial prediction with prob = 0
  if(length(unique(labels)) == 1)
  {
    if(sign(labels[1]) == -1)
      labels <- c(labels,1)
    else
      labels <- c(labels,-1)
    prob <- c(prob,0)
  }

  pred <- prediction(prob, labels)
  unlist(performance(pred, "auc")@y.values)
}
