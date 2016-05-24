#' #############################################################################################
#'  Cross validation for FrSVM, an R algorithm, which integrates protein-protein 
#'		 		   interaction network information into gene selection for microarry classification
#'
#' @Parameters---- 
#' @param x: gene expression data
#' @param y: class labels
#' @param d: damping factor for GeneRank, defaults value is 0.5
#' @param Gsub: Adjacency matrix of Protein-protein intersction network
#' @param folds: # of -folds cross validation (CV)
#' @param repeats: number of CV repeat times
#' @param parallel: paralle computing or not
#' @param cores: cores used in parallel computing
#' @param DEBUG: show more results or not
#' @param Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @seed  seed: for random sampling.
#' @param top.uper: the uper bound of top ranked genes
#' @param top.lower: the lower bound of top ranked genes
#'
#' @Returned results---- 
#' \item{auc}{AUC values of each test fold in CV}
#' \item{labels}{original class labels}
#' \item{fits}{SVM model for each training fold}
#' \item{feat}{Selected features in each training folds}

#' @references Yupeng Cun, Holger Frohlich (2012) Integrating Prior Knowledge Into Prognostic Biomarker Discovery Based on Network Structure.arXiv:1212.3214 
#' @references Winter C, Kristiansen G, Kersting S, Roy J, Aust D, et al. (2012) Google Goes Cancer: Improving Outcome Prediction for Cancer Patients by Network-Based Ranking of Marker Genes. PLoS Comput Biol 8(5): e1002511. doi:10.1371/journal.pcbi.1002511
#' @export
#' library(netClass)
#' 
#' data(expr)
#' data(ad.matrix)
#' x <- expr$genes
#' y <- expr$y
#' 
#' 
# r.frsvm <- cv.frsvm(x=x, y, folds=5,Gsub=ad.matrix, repeats=3, parallel = FALSE, cores = 2, DEBUG=TRUE,d=0.85,top.uper=10,top.lower=50,seed=1234,Cs=10^c(-3:3))


cv.frsvm <- function(x, y, folds=10,Gsub=matrix(1,100,100), repeats=5, parallel = FALSE, cores = 2, DEBUG=FALSE,d=0.85, top.uper=10,top.lower=50,seed=1234, Cs=10^c(-3:3))
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
  
	
	if (folds < 2) stop("folds should be greater than or equal to 2.\n")
	if (folds > n) stop("folds should be less than or equal to the number of observations.\n")

	cuts  <- cv.repeats <- list()  	  
	op =	top.uper 
	aa =    top.lower 
	
	int= intersect(colnames(x),colnames(Gsub))
	x=x[,int]
	Gsub=Gsub[int,int]
			
	set.seed(1234)	
	for(r in 1:repeats)
	{
		perm = sample(1:n)		
		#perm <- sample(1:n) #Sampling a random integer between 1:n
		repeat.models <- NULL 
    
		for(k in 1:folds) #randomly divide the training set in to 10 folds
		{
			tst <- perm[seq(k, n, by=folds)]  #      
			trn <- setdiff(1:n, tst)            
			cuts[[k]] <- list(trn=trn, tst=tst)    
		}    	
    
		if(DEBUG) cat('Starting classification of repeat:',r,'\n')
    
		if(parallel)	repeat.models <- mclapply(1:folds, classify.frsvm, cuts=cuts, x=x, y=y,cv.repeat=r, DEBUG=DEBUG,Gsub=Gsub,d=d,op=op,aa=aa,Cs=Cs)
		else		repeat.models <-   lapply(1:folds, classify.frsvm, cuts=cuts, x=x, y=y, cv.repeat=r, DEBUG=DEBUG, Gsub=Gsub,d=d,op=op,aa=aa,Cs=Cs)
	
		
    
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
  	
  	feat <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$feat))  
  	names(feat) <- paste("Repeat",1:repeats,sep="")  
  	feat <- lapply(feat, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })
  	
  	res <- list(feat=feat, auc=auc,fits=fits, labels=y)  
  	class(res) <- 'netClassResult' 
  	
  	return(res)  
	#return ( cv.repeats)
}



# Training and predicting using FrSVM classification methods
#-----------------------------  
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels.
# DEBUG: show debugging information in screen or not.
# d: damping factor for GeneRank, defaults value is 0.5
# op: the uper bound of top ranked genes
# aa: the lower bound of top ranked genes
# Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
# Gsub: an adjacency matrix that represents the underlying biological network.
#-----------------------------  
# Return a list with the results of the training model and predcit AUC   


classify.frsvm <- function(fold, cuts, x, y, cv.repeat, DEBUG=DEBUG,Gsub=Gsub,d=d,op=op,aa=aa,Cs=Cs)
{
	gc() 	
	if(DEBUG) cat('starting Fold:',fold,'\n')	
  	#ad.list<- as.adjacencyList(Gsub)
	## get training and test indices
	trn <- cuts[[fold]]$trn
	tst <- cuts[[fold]]$tst		

	
	train	<- train.frsvm(x=x[trn,], y=y[trn], DEBUG=DEBUG,Gsub=Gsub,d=d,op=op,aa=aa,Cs=Cs)	    	  	  
	test	<- predictFrsvm(train=train, x=x[tst,], y=y[tst], DEBUG=DEBUG)
	
	auc		<- test$auc
	feat	<- train$feat						
	
	if(DEBUG) cat("=> the best AUC   is ", auc, "\t Best features length:  ", length(feat),"\n")					
	if(DEBUG) cat('Finished fold:',fold,'\n\n')

	gc()	
	res=list(fold=fold, model=train, auc=auc,feat= feat)
		
	return(res)
}

train.frsvm <- function(x=x, y=y, DEBUG=FALSE,Gsub=Gsub,d=0.85,op=10,aa=50,Cs=10^(-3:3))
{		    	  	  
	fits  <- list()
  	best.boundL <- list()
  	featL  <- list()
 
	if(DEBUG) cat("Geting Gene Ranking \n")			
	ranks = getGeneRanking(x=x, y=y, Gsub=Gsub, d=d)			
	
	ranksI=sort(ranks,decreasing=T)
	#nodesI=sort(topNodes,decreasing=T)
			
	i=1
	nn=aa		
	while(nn > op)
	{ 		  				  		
			feat.rank = names(ranksI[1:nn])	  	
	  		feat = feat.rank
	  		featL[[i]] = feat			 
		  	#print(length(feat))	
	  		fit <- svm.fit(x=x[,feat], y=y, Cs=Cs, scale="scale", DEBUG=FALSE) 	  		  
	  		fits[[i]]=fit	  		  		
		    best.boundL[[i]] <-  fits[[i]]$error.bound  
			nn=nn-3
			i=i+1
	}	
		
	if(DEBUG)cat("the opitimal steps: ", i-1, "\n")
		
	best.boundLs= unlist(best.boundL)
	best.index = which(best.boundLs==min(best.boundLs))		
	n=length(best.index)						
			
	trained   =	fits[[best.index[n]]]		
	feat = featL[[best.index[n]]]	
	
	trained$graph = graph.adjacency(Gsub[feat,feat], mode="undirected")	
	
	res=list(trained=trained, feat= feat)
		
	return(res)
}

predictFrsvm <- function(train=train, x=x, y=y, DEBUG=FALSE)
{
	feat = train$feat	
	xts= x[,feat]				
	label 	<- sign(as.numeric(y) - 1.5) # for computing the AUC
	test    <- svm.predict(fit=train$trained, newdata=xts, type="response")	
	## calculate the AUC					
	auc	<- calc.auc(test, label)							
	#if(DEBUG) cat("=> the best AUC   is ", auc, "\t Best features length:  ", length(feat),"\n")				
	res=list(auc=auc,pred=test)	
	return(res)
}

## geting gene ranking of differencial expression of <xi,y> based on related network
getGeneRanking = function(x=x, y=y,Gsub=Gsub, d=d)
{
	
    int= intersect(colnames(x),colnames(Gsub))
    x=x[,int]
    Gsub=Gsub[int,int]     
    x = scale(x)
    #calculate the t-scorce of each probe.
    xtt = matrix(0, ncol(x))
    yy=sign(as.numeric(y)-1.5)

    for(i in 1: ncol(x))
    	xtt[i]= abs(as.numeric(t.test(x[,i], yy,paired=TRUE)$statistic))
    
    names(xtt)= colnames(x)	
    exprs = xtt[,1]      	
    names(exprs) = colnames(x)
    ranks   <- pGeneRANK(W=Gsub, ex=exprs, d=d)   

    names(ranks) <- colnames(Gsub)  	
    return(ranks)
}



## GeneRank from pathClas 
pGeneRANK <- function(W,ex,d, max.degree=Inf)
{
  #require(Matrix)
 #W=Matrix(W)
 
  ex = abs(ex)

  ## normalize expression values
  norm_ex = ex/max(ex)

  ## try sparse Matrices in R => later
  ##w = sparse(W)
  dimW = dim(W)[1]
  if(dim(W)[2]!=dimW) stop("W must be a square matrix.")
   
  ## get the in-dgree for every node
  ## from KEGG we get a directed graph
  ## thus, the column sums correspond to
  ## to the in-degree of a particular gene
 
  degrees = pmin(max.degree, pmax(1,colSums(W), na.rm=T))

  ## A = Identity Matrix with dimensions
  ## same as the adjacency matrix

  A=Matrix(0, nrow = dimW, ncol = dimW)
  diag(A) = 1
  
  ## produce a matrix with the degrees on
  ## the diagonal

  D1=Matrix(0, nrow = dimW, ncol = dimW)
  diag(D1) = 1.0/degrees

  ## divide the in-degrees of the gene
  ## by the overall in-degree(colSum) of the gene
  ## => kind of normalizing the in-degrees

  A = A - d*(Matrix(t(W)) %*% D1)

  ## here, we give 1-d 'for free'
  b  = (1-d) * norm_ex

  ## we want to solve:
  ## (I - d W^t D^-1)r = (1-d)ex which is the Jacobi of the PageRank
  ## where A = (I - d W^t D^-1)
  ## and   b = (1-d)ex
  ## therefore:  Ar = b
  r = as.numeric(solve(A,b))
  return(r)
  
}
