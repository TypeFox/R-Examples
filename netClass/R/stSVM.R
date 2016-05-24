#' Cross validation for smoothed t-statistic to select significant top ranked differential expressed genes
#'
#'
#' @param x: a p x n matrix of expression measurements with p samples and n genes.
#' @param y: a factor of length p comprising the class labels.
#' @param DEBUG: show debugging information in screen or not.
#' @param Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @param Gsub: an adjacency matrix that represents the underlying biological network.
#' @param pt.pvalue: cut off p value of permutation test
#' @param aa:  permutation test steps
#' @param op: optimal on top op% ranked genes
#' @param a: constant value of random walk kernel
#' @param p: random walk step(s) of random walk kernel
#' @param allF: Using all features (TRUE) or only these genes mapped to prior information (FALSE). 
#' @seed  seed: seed for random sampling.
#' @return a list with the results of the cross-validation.  \item{feat}{the selected features}  \item{fit}{the fitted SVM model} \item{auc}{the prediction auc of the fitted SVM model} \item{lables}{the test labels}  
#' @references Y.Cun and H. Frohlcih, (2013), Data Integration via Network Smoothed T-Statistics, submitted    
#' @export
#' @note The optimal number of features is found by using the span estimate. See Chapelle, O., Vapnik, V., Bousquet, O., and Mukherjee, S. (2002). Choosing multiple parameters for support vector machines. \emph{Machine Learning}, 46(1), 131-159.
#'
#' @author  Yupeng Cun \email{yupeng.cun@gmail.com}
#' @examples
#' 
#' data(expr)
#' x <- expr$genes
#' y <- expr$y
#'
#' library(netClass)
#' 
#' r.stsvm <- cv.stsvm(x=x,x.mi=NULL,y=y, folds=5,Gsub=ad.matrix,  repeats=3, parallel = TRUE, cores=4,DEBUG=TRUE,pt.pvalue=0.05,op=0.9,aa=1000,a=1,p=2,allF=TRUE,Cs=10^(-4:4))
#' mean(r.stsvm$auc)


cv.stsvm <- function(x=x, x.mi=NULL,y=y, folds=5,Gsub=matrix(1,100,100),op.method=c("pt","spb"), repeats=3, parallel = FALSE, cores=2,DEBUG=TRUE, pt.pvalue=0.05,op=0.85, aa=1000,a=1,p=2,allF=TRUE,seed=1234,Cs=10^c(-3:3))
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
	if(is.null(x.mi))
	{
		x= cbind(x,x.mi)
	}
	
  	dk.tf = pt.pvalue
  	pred.type="class"
  	#aa=pt.step
  	mapping=NULL
	ex.sum=x
  	
	if(allF)
	{
		
		int.all = colnames(x)
		cat("All ",length(int.all),"features.\n")
		ad.all <- matrix(0, ncol=length(int.all), nrow=length(int.all))
		colnames(ad.all) <- colnames(x)
		rownames(ad.all) <- colnames(x)
		inta = intersect(int.all, colnames(Gsub))
		ad.all[inta,inta]= Gsub[inta,inta]		
		Gsub <- ad.all
	}	
  
	
	int = intersect(colnames(Gsub),colnames(x))
	x=x[,int]
	Gsub=Gsub[int,int]
	
	#Calculating random walk kernel matrix of network.
	dk <- calc.diffusionKernelp(L=Gsub, is.adjacency=TRUE,p=p,a=a)
	
	#if(fs.method=="TRG+SVM")dk <- calc.diffusionKernelp(L=Gsub[int,int], is.adjacency=TRUE,p=p,a=a)
	
	cuts  <- cv.repeats <- list()  
				
	set.seed(1234)	
	for(r in 1:repeats)
	{		
		###random sampling 
		stratified=FALSE
		if(stratified){
			perm0 = sample(which(y == levels(y)[1]))
			perm1 = sample(which(y == levels(y)[2]))
			perm = c(perm0, perm1)
		}
		else
			perm = sample(1:n)		
		#perm <- sample(1:n) #Sampling a random integer between 1:n
		repeat.models <- NULL 
    
		for(k in 1:folds) #randomly divide the training set in to 10 folds
		{
			tst <- perm[seq(k, n, by=folds)]  #      
			trn <- setdiff(1:n, tst)            
			cuts[[k]] <- list(trn=trn, tst=tst)    
		}    
		
		#pb <- txtProgressBar(min = 0, max = folds, style = 3)   #txtProgressBar id funtion in packages "utils", Text progress bar in the R console
    
		if(DEBUG) cat('Starting classification of repeat:',r,'\n')
    
		if(parallel)	repeat.models <- mclapply(1:folds, classify.stsvm, cuts=cuts, ex.sum=ex.sum, x=x, y=y,cv.repeat=r, DEBUG=DEBUG, Gsub=Gsub,op.method=op.method,op=op,aa=aa,dk=dk,dk.tf=dk.tf,p=p,a=a,seed=seed,Cs=Cs)
		else		repeat.models <-   lapply(1:folds, classify.stsvm, cuts=cuts, ex.sum=ex.sum,x=x, y=y,cv.repeat=r, DEBUG=DEBUG, Gsub=Gsub,op.method=op.method,op=op,aa=aa,dk=dk, dk.tf=dk.tf,p=p,a=a,seed=seed,Cs=Cs)
	
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
  	
  	feat <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$feat))  
  	names(feat) <- paste("Repeat",1:repeats,sep="")  
  	feat <- lapply(feat, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })
  	
  	res <- list(feat=feat, auc=auc,fits=fits, labels=y)  
  	class(res) <- 'netClassResult' 
  	#res<-cv.repeats
  	return(res)  
	#return ( cv.repeats)
}

# Training and predicting using stSVM classification methods
#-----------------------------  
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels.
# DEBUG: show debugging information in screen or not.
# Cs: soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
# Gsub: an adjacency matrix that represents the underlying biological network.
# pt.pvalue: cut off p value of permutation test
# aa:  permutation test steps
# op: optimal on top op% ranked genes
#-----------------------------  
# Return a list with the results of the training model and predcit AUC   


classify.stsvm <- function(fold, cuts,ex.sum, x, p,a, y, cv.repeat, DEBUG=DEBUG,Gsub=Gsub,op.method=op.method, op=op,aa=aa,dk=dk,dk.tf=dk.tf,seed=seed,Cs=Cs)
{
	gc() 		
	
	if(DEBUG) cat("stSVM ==>> starting Fold:",fold,"\n")	
	
  
	## get training and test indices
	trn <- cuts[[fold]]$trn
	tst <- cuts[[fold]]$tst		
	#print(tst)
	label 	<- sign(as.numeric(y[tst]) - 1.5) # because y is a factor			
	
	train	<- train.stsvm(x=x[trn,], y=y[trn], DEBUG=DEBUG,Gsub=Gsub, op.method=op.method, op=op,aa=aa,dk=dk, dk.tf=dk.tf,seed = seed,Cs=Cs)	    	  	  
	test	<- predictStsvm(train=train, x=x[tst,], y=y[tst],DEBUG=DEBUG)
	
	auc		<- test$auc
	feat	<- train$feat						
	
	if(DEBUG) cat("=> the best AUC   is ", auc, "\t Best features length:  ", length(feat),"\n")					
	if(DEBUG) cat('Finished fold:',fold,'\n\n')

	gc()	
	res=list(fold=fold, model=train, auc=auc,feat= feat)
		
	return(res)
}


train.stsvm <- function(x=x, y=y, DEBUG=FALSE,Gsub=Gsub, op.method="sp", op=10,aa=100,dk=dk, dk.tf=0.05,seed = 1234,Cs=10^(-3:3), EN2SY=NULL)
{
	# ranked by smoothed t-statistics  	
	sca=TRUE			
	ranks = getGraphRank(x=x, y=y, Gsub=dk,sca=sca)							
	
	
	if(op.method=="pt")
	{
		topGenes = names(ranks[which(ranks > quantile(ranks,op))])				
		topRanks =ranks[topGenes]
		#if(DEBUG) cat("range of topRanks:",range(topRanks),"\n=======================\n")
		
		ytrn = y
		xtrn =x[,topGenes]			

		R = aa
		set.seed(seed)
		p_ranks=matrix(0,ncol=length(topRanks))[1,]
		names(p_ranks)=topGenes
		
		if(DEBUG) cat("\nStart permutation test for significant topRanks. \n")	
		
		for(i in 1 :R )
		{							
			ys = sample(ytrn)															
			rankt <-   getGraphRank(x=xtrn, y=ys, Gsub=dk, sca=sca)			
			for(j in topGenes) if(rankt[j] >= topRanks[j])p_ranks[j] = p_ranks[j]+1					
		}	
		cat("head(p_ranks)   ", head(p_ranks),"\n")
		p_ranks = p_ranks/R	
		#p_ranks = p.adjust(p_ranks, method="BH")
		names(p_ranks)=names(topRanks)
	
		sigRank <- names(p_ranks[which(p_ranks < dk.tf)])		
		if(length(sigRank)==0)sigRank=names(sort(topRanks))[1:10]									
		feat  = sigRank	
		train 	 <- svm.fit(x=x[,feat], y=y, Cs=Cs, scale="scale", DEBUG=FALSE) 					
	}
	
	if(op.method=="sp")
	{
		fits  <- list()
	  	best.boundL <- list()
	  	featL  <- list()
	 
		if(DEBUG) cat("Geting Gene Ranking \n")			
		#ranks = getGeneRanking(x=x, y=y, Gsub=Gsub, d=d)			
	
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
			
		train   =	fits[[best.index[n]]]		
		feat = featL[[best.index[n]]]	
	}	

	if ( !is.null(EN2SY))
	{		
		syb = EN2SY[feat,2]
		sGsub=Gsub[feat,feat]
		colnames(sGsub) = syb
		rownames(sGsub) = syb
		train$graph = graph.adjacency(sGsub, mode="undirected")		
	}
	else
	{
		train$graph = graph.adjacency(Gsub[feat,feat], mode="undirected")	
	}
	train$feat=feat	
	
	res=list(trained=train, feat= feat)		
	return(res)
}



predictStsvm <- function(train=train, x=x, y=y,DEBUG=DEBUG)
{	
	label 	<- sign(as.numeric(y) - 1.5) # because y is a factor			
	
	feat  = train$feat				
	test    <- svm.predict(fit=train$trained, newdata=x, type="response")		
	auc	<- calc.auc(test, label)				
	res=list(auc=auc,pred=test)	
	return(res)
}


## 
# Computing the Random Walk Kernel matrix of network
#-----------------------------  
# L: an adjacency matrix that represents the underlying biological network.
# a: constant value of random walk kernel
# p: #(p) random walk step(s) of random walk kernel
#-----------------------------  

calc.diffusionKernelp = function(L, is.adjacency=TRUE, p=3,a=2)
{
  if(missing(L)) stop('You have to provide either a laplace or a adjacency matrix to calculate the diffusion kernel.')
  print("thresh")
  
  if(is.adjacency)
  {    
    dnames <- dimnames(L)
    L = graph.adjacency(L, mode='undirected', diag=FALSE)
    L = graph.laplacian(L, normalized=TRUE)
    dimnames(L) = dnames
  }
  
  n=ncol(L)
  I=diag(n)      
      
  if( p==1) R = a*I-L
    
  else
  {
    R=a*I -L
    for(ii in seq(from=2,to=p, by=1))
    {
      R = R %*% (a*I -L)
      ii = ii+1
    }
  }
  
  #KernelDiag <- sqrt(diag(R) + 1e-10)
  #R.norm 	<-  R/(KernelDiag %*% t(KernelDiag))
  #R			<-  R.norm  
  
  colnames(R) =colnames(L)
  rownames(R) =rownames(L)
  R    
}

##### random walk kernel matrix smoothing t-statistic
#------------------------------------
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels. 
# Gsub: Random Walk Kernel matrix of network
# -----------------------------------------
getGraphRank = function(x=x, y=y, Gsub=Gsub,sca=TRUE)
{    
    int= intersect(colnames(x),colnames(Gsub))  
    x=x[,int]
    Gsub=Gsub[int,int]   
    #print(length(int)) 
    if(sca==TRUE) x = scale(x)	
    #calculate the t-scorce of each probe.
    xtt = rep(0, ncol(x))
    yy=sign(as.numeric(y)-1.5)
    yy[yy==-1]=0

    for(i in 1: ncol(x))  xtt[i]= t.test(x[,i], yy,paired=T)$statistic    	   
    
    names(xtt)= colnames(x)	
    exprs = abs(xtt)
    exprs = exprs/ max(exprs)       	
    ranks   <- t(exprs[int]) %*% Gsub[int,int]
    r=ranks[1,]
    r= r/max(r)
    #print(dim(ranks))
   # print(dim(Gsub))
    names(r) = colnames(x)	 	   
   # print("finised GrapRank")
    return(r)
}				


