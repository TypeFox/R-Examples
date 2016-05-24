#' Cross validation for hub nodes classification
#'
#'
#' @param x: a p x n matrix of expression measurements with p samples and n genes.
#' @param y: a factor of length p comprising the class labels.
#' @param DEBUG: show debugging information in screen or not.
#' @param scale: a character vector defining if the data should be centered and/or scaled.
#' Possible values: are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
#' @param Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @param Gsub: an adjacency matrix that represents the underlying biological network.
#' @param Gs: Undirected of graph with adjacency matrix Gsub.
#' @param node.ct: cut off value for select highly quantile nodes in a nwtwork. Defaults to \code{0.98)}.
#' @seed  seed: for random sampling.
#' @return a list with the results of the cross-validation.  \item{feat}{the selected features}  \item{fit}{the fitted SVM model} \item{auc}{the prediction auc of the fitted SVM model} \item{lables}{the test labels}  
#' @references Taylor et al.(2009)Dynamic modularity in protein interaction networks predicts breast cancer outcome, Nat. Biotech.: doi: 10.1038/nbt.1522
#' @export
#' @note The optimal number of features is found by using the span estimate. See Chapelle, O., Vapnik, V., Bousquet, O., and Mukherjee, S. (2002). Choosing multiple parameters for support vector machines. \emph{Machine Learning}, 46(1), 131-159.
#'
#' @author  Yupeng Cun \email{yupeng.cun@gmail.com}
#' @examples
#' library(netClass)
#' data(expr)
#' data(ad.matrix)
#' Gs = graph.adjacency(ad.matrix[1:1000,1:1000], mode="undirected")
#' x <- expr$genes
#' y <- expr$y
#' r.hubC <- cv.hubc(x=x, y=y, folds=10, repeats=3, parallel = TRUE, cores = 4, DEBUG=TRUE,nperm=1000,Gsub=ad.matrix,Gs=Gs2,node.ct=0.95,Cs=10^(-3:3))

cv.hubc <- function(x, y, folds=10, repeats=5, parallel = TRUE, cores = NULL, DEBUG=TRUE, nperm=500,  node.ct=0.98, Gsub=matrix(1,100,100), Gs=Gs,seed=1234,Cs=10^c(-3:3))
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

	n     <- length(y)
	folds <- trunc(folds)
	if (folds < 2) stop("folds =< 2.\n")
	if (folds > n) stop("folds > the number of observations.\n")
		
	
	### select the signaficant hubs
		
	int = intersect(colnames(Gsub), graph::nodes(Gs))
	int = intersect(int, colnames(x))
	x=x[,int]	
	Gsub = Gsub[int,int]
	
	hubs <- which(colSums(Gsub) > quantile(colSums(Gsub), node.ct))
	degM <- Gsub[hubs,hubs]
	#hubs <- colnames(degM)
	#print(hubs)
	
	int = intersect(colnames(Gsub), colnames(degM))
	hubs = intersect(hubs, colnames(x))
	#x=x[,int]	
	#Gsub = Gsub[int,int]
	gHub = subGraph(hubs,Gs)
	#print(hubs)
	
	nHub = length(hubs)
	  
	
	
	cuts  <- cv.repeats <- list()  
	set.seed(seed)
	for(r in 1:repeats)
	{
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
    
		if(parallel)	repeat.models <- mclapply(1:folds, r=r, classify.hubc, cuts=cuts, x=x, y=y, cv.repeat=r, Gsub=Gsub,DEBUG=DEBUG,gHub=gHub, hubs=hubs, nperm=nperm,node.ct=node.ct,Cs=Cs)
		else		repeat.models <-   lapply(1:folds,r=r, classify.hubc, cuts=cuts, x=x, y=y, cv.repeat=r, Gsub=Gsub,DEBUG=DEBUG,gHub=gHub, hubs=hubs, nperm=nperm,node.ct=node.ct,Cs=Cs)
	
		#close(pb)
    
		if(length(repeat.models) != folds)
		{      
			geterrmessage()      
			stop("One or more processes did not return. May be due to lack of memory.\n")
		}
	
		if(DEBUG)	cat('All models of repeat:',r,'have been trained.\n')
    
		cv.repeats[[r]] <- repeat.models
  
  	} 
  	
  	feat <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$feat))  
  	names(feat) <- paste("Repeat",1:repeats,sep="")  
  	feat <- lapply(feat, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })
  		

  	auc <- sapply(cv.repeats, function(cv.repeat) sapply(cv.repeat, function(model) model$auc))  
  	colnames(auc) <- paste("Repeat",1:repeats,sep="")  
  	rownames(auc) <- paste("Fold",1:folds,sep="")  
  	
  	fits <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$model))  
  	names(fits) <- paste("Repeat",1:repeats,sep="")  
  	fits <- lapply(fits, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })  
  	

  	res <- list(feat=feat, auc=auc, fits=fits, labels=y)  
  	class(res) <- 'netClassResult' 
  	return(res)  
	#return ( cv.repeats)
}


# Training and predicting using hub nodes classification methods
#-----------------------------  
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels.
# DEBUG: show debugging information in screen or not.
# hubs: hub nodes
# nperm: the repeat number of permutation test 
#-----------------------------  
# Return a list with the results of the training model and predcit AUC   


classify.hubc <- function(fold, r, cuts, x, y, cv.repeat,Gsub=Gsub, DEBUG=DEBUG, gHub=gHub, hubs=hubs,nperm=nperm,node.ct=node.ct,Cs=Cs)
{
	gc()  
	
	if(DEBUG) cat("\nstarting reapet: #",r," # of Fold: #",fold,'\n')
  
	## get training and test indices
	trn <- cuts[[fold]]$trn
	tst <- cuts[[fold]]$tst

	train	<- train.hubc(x=x[trn,], y=y[trn], DEBUG=DEBUG, Gsub=Gsub, gHub=gHub, hubs=hubs,nperm=nperm,node.ct=node.ct,Cs=Cs)	    	  	  
	test	<- predictHubc(train=train, x=x[tst,], y=y[tst], DEBUG=DEBUG)
	
	auc		<- test$auc
	feat	<- train$feat						
	
	if(DEBUG) cat("=> the best AUC   is ", auc, "\t Best features length:  ", length(feat),"\n")					
	if(DEBUG) cat('Finished fold:',fold,'\n\n')

	gc()	
	res=list(fold=fold, model=train, auc=auc,feat= feat)
		
	return(res)
}

train.hubc <- function(x=x, y=y, DEBUG=FALSE, Gsub=Gsub, gHub=gHub, hubs=hubs,nperm=500,node.ct=0.95,Cs=10^(-3:3))
{
	
	if(DEBUG)cat("\nlength(hubs):",length(hubs),"\n")
	p_hub <- pOfHubs(x=x,y=y, gHub=gHub, hubs=hubs, nperm=nperm)	
	index = p_hub$hub[which(p_hub$pVal < 0.5)]		
	#if(length(index) <=1)
	#{
	#  tt=p_hub$pVal
	#  names(tt)=p_hub$hub
	#  tnodes = p_hub$hub[which( tt > quantile(tt, node.ct))]
	#  index=  p_hub$hub[which(p_hub$hub %in% tnodes)]
	#}
	
	if(length(index)>=1) feat = index
	if(length(index)<1)	feat=hubs
	trained <- fit.pe(x=x[,feat],y= y, DEBUG=DEBUG, scale="center", Cs=Cs)  
	feat=colnames(x[,feat])
	## save the test indices
	##trained$graph = graph.adjacency(Gsub[feat,feat], mode="undirected")	
		
	res=list(trained=trained, feat= feat)		
	return(res)
}

predictHubc <- function(train=train, x=x, y=y,DEBUG=FALSE)
{
	feat = train$feat				
	label 	<- sign(as.numeric(y) - 1.5) # for computing the AUC  
	test    <- svm.predict(fit=train$trained$fit, newdata=x,type="response")
	
	auc   <- calc.auc(test, label)
	
	if(DEBUG) cat("=> the best AUC   is ", auc, "\t Best features length:  ", length(feat),"\n")				
	res=list(auc=auc,pred=test)	
	return(res)
}

##################################
###
### Calculated the significant Hubs
###

XIs <- function(xm=exprs, hub=hub, gHub=gHub)
{
	inters = unlist(adj(gHub, hub))	
	
	nColx= ncol(xm)
	inters = inters[which(inters<=nColx)]
	n=length(inters)
	xis= matrix(0, nrow(xm), n)
	
	for (i in 1:n)
	{			
		xis[,i]=xm[,inters[i]]  #cat("\n :", i, "\t",inters[i],"\n" )		
	}		
	return (xis)
}
#XIs(xm=xm, hub=hub, gHub=gHub)

XH <- function(xm=exprs, hub=hub, gHub=gHub)
{
	xh= matrix(0, nrow(xm))
	
	xh=xm[,hub]
		
	return (xh)
}

aveHubDiff <- function(xA=xA, xD=xD, gHub=gHub, hub=hub)
{	
	inters = unlist(adj(gHub, hub))	
	nColx= ncol(xA)+ncol(xD)
	inters = inters[which(inters<=nColx)]
	n=length(inters)		
	
	AverageHubDiff=0.0
	
	if(n!=0)
	{
		aveHubDiff= matrix(0, n)	
		hA = XH(xm=xA,hub=hub,gHub=gHub)
		hD = XH(xm=xD,hub=hub,gHub=gHub)
		hM = (mean(hA)+ mean(hD))/2.0
	
		iAs = XIs(xm=xA,hub=hub,gHub=gHub)
		iDs = XIs(xm=xD,hub=hub,gHub=gHub)
		for(i in 1:n)
		{
			iA = iAs[,i]
			iD = iDs[,i]
			iM = (mean(iA)+ mean(iD))/2.0
			
			na = nrow(xA)
			nd = nrow(xD)
	
			sHIa = sum((iA-iM)*(hA-hM))	
			sHId = sum((iD-iM)*(hD-hM))
	
			varIa = sqrt( sum((iA-iM)*(iA-iM)) / (na-1) )
			varHa = sqrt( sum((hA-hM)*(hA-hM)) / (na-1) )
			varId = sqrt( sum((iD-iM)*(iD-iM)) / (nd-1) )
			varHd = sqrt( sum((hD-hM)*(hD-hM)) / (nd-1) )
	
			deltaRADi= sHIa/((na-1)*varIa*varHa) - sHId/((nd-1)*varId*varHd)
		
			aveHubDiff[i] = abs(deltaRADi)
			#cat("aveHubDiff[",i,"]:	", aveHubDiff[i],"\n")
		}	
		AverageHubDiff= sum(aveHubDiff)/n
		#cat("AverageHubDiff", AverageHubDiff,"\n")
	}
	return (AverageHubDiff)
}#ahd <- aveHubDiff(xA=xA, xD=xD, gHub=gHub, hub=hub)

#Random Permutation Test
randomPermutation.test <- function(x=x,y=y, gHub=gHub, hub=hub, nperm=nperm) 
{
	xA=x[y==1,]
	xD=x[y==-1,]
	avhd.orig = aveHubDiff(xA=xA, xD=xD, gHub=gHub, hub=hub)
	p = 0
	set.seed(1234)
	for(i in 1:nperm)
	{
		ys = sample(y)
		avhd = aveHubDiff(x[ys==1,], x[ys==-1,],gHub=gHub, hub=hub)
		#cat("\n\t aveHubDiff:", avhd, "\n")
		if(avhd >= avhd.orig)
			p = p + 1		
	}
	p=p / nperm

	return (p)
}
#randomPermutation.test(x=xm,y=y, gHub=gHub, hub=hub, nperm=nperm)


pOfHubs<- function(x=x,y=y, gHub=gHub, hubs=hubs, nperm=nperm) 
{
	pHubs = matrix(0, length(hubs))	
	
	for(i in 1: length(hubs))
	{
		hub = hubs[i]
		cat("\n",hub," ")
		pHubs[i] = randomPermutation.test(x=x,y=y, gHub=gHub, hub=hub, nperm=nperm)
		cat("p[",i,"] vaule of hub [",hub,"] after ",nperm, " repeats: ", pHubs[i])	
	}
	#pHubs=p.adjust(pHubs, method="BH")
	p_hub= list(pVal=pHubs, hub=hubs )
	return(p_hub)	
}



######################
# kernel bound SVM 
###

fit.pe <- function(x, y, DEBUG=FALSE, scale=c('center', 'scale'), Cs=10^c(-3:3))
{ 
	best.bound = Inf
  	feat = colnames(x)	
	fit = svm.fit(x=x, y=y, Cs=Cs, scale="center", DEBUG=TRUE)      
          
	if(DEBUG) cat('Best Model is: Spanbound=',fit$error.bound,', C=',fit$C,',', length(fit$features),'features.\n')
  
	result <- list(features=feat, error.bound = fit$best.bound, fit=fit)
  	#class(result) <- 'sg'
	
	return(result)
}


