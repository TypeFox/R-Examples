#' Cross validation for Pathway Activities Classification(PAC)
#'
#' Implementation of the Pathway Activities Classification  aepSVM) by CROG algorithm.
#'
#' @param x: a p x n matrix of expression measurements with p samples and n genes.
#' @param y: a factor of length p comprising the class labels.
#' @param DEBUG: show debugging information in screen or not.
#' @param Gsub: an adjacency matrix that represents the underlying biological network.
#' @seed  seed: for random sampling.
#' @return a list with the results of the cross-validation.  \item{feat}{the selected features}  \item{fit}{the fitted SVM model} \item{auc}{the prediction auc of the fitted SVM model} \item{lables}{the test labels}  
#' @references Lee E, Chuang H-Y, Kim J-W, Ideker T, Lee D (2008) Inferring Pathway Activity toward Precise Disease Classification. PLoS Comput Biol 4(11): e1000217. doi:10.1371/journal.pcbi.1000217
#' @export
#' 
#' @author  Yupeng Cun \email{yupeng.cun@gmail.com}

#' @examples
#' library(netClass)
#' 
#' data(expr)
#' data(ad.matrix)
#' x <- expr$genes
#' y <- expr$y
#'
#' r.pac <- cv.pac(x=x, y=y, folds=10, repeats=10, parallel = FALSE, cores = 4, DEBUG=TRUE, Gsub=ad.matrix)

cv.pac <- function(x=x, y=y, folds=10, repeats=5, parallel = TRUE, cores = NULL, DEBUG=TRUE, Gsub=matrix(1,100,100), seed=1234)
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
	int= intersect(colnames(Gsub), colnames(x))
	Gsub = Gsub[int,int]	
	x=x[,int]
	
	if (folds < 2) stop("folds should >= 2.\n")
	if (folds > n) stop("folds should > the number of observations.\n")
		
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
    
		if(parallel)	repeat.models <- mclapply(1:folds, classify.pac, Gsub=Gsub, cuts=cuts,x=x, y=y, cv.repeat=r, int=int,DEBUG=DEBUG)
		else		repeat.models <-   lapply(1:folds, classify.pac, Gsub=Gsub,cuts=cuts, x=x, y=y, cv.repeat=r,int=int, DEBUG=DEBUG)
	
		#close(pb)
    
		if(length(repeat.models) != folds)
		{      
			geterrmessage()      
			stop("One or more processes lost. May be due to lack of memory.\n")
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
  	
  	
  	res <- list(feat=feat, auc=auc, fits=fits, labels=y)  
  	class(res) <- 'netClassResult' 
  	return(res)  
	#return ( cv.repeats)
}



# Training and predicting using PAC classification methods
#-----------------------------  
# x: a p x n matrix of expression measurements with p samples and n genes.
# y: a factor of length p comprising the class labels.
# DEBUG: show debugging information in screen or not.
#-----------------------------  
# Return a list with the results of the training model and predcit AUC   


classify.pac <- function(fold, cuts, x, y, cv.repeat,Gsub=Gsub,int, DEBUG=FALSE)
{
	gc()  	
	if(DEBUG) cat('starting Fold:',fold,'\n')
  
	## get training and test indices
	trn <- cuts[[fold]]$trn
	tst <- cuts[[fold]]$tst	
	
	train	<- train.pac(x=x[trn,], y=y[trn], int=int, DEBUG=DEBUG,Gsub=Gsub)	    	  	  
	test	<- predictPac(train=train,x=x[tst,], y=y[tst],int=int, DEBUG=DEBUG)
	
	auc		<- test$auc
	feat	<- train$crog						
	
	if(DEBUG) 
	{
		cat("=> the best AUC   is ", auc, "\t Best features length:  ", length(feat),"\n")					
		cat('Finished fold:',fold,'\n\n')
	}

	gc()	
	res=list(fold=fold, model=train, auc=auc,feat= feat)
		
	return(res)
}


train.pac <- function(x=x, y=y, int=int, DEBUG=FALSE,Gsub=Gsub)
{	
	yy=sign(as.numeric(y))
	yy[yy==-1]=0
	
	apTrain<-probeset2pathwayTrain(x=x,y=y, int=int)		
	mydata = data.frame(ap=apTrain$ap,yy=yy)		
	trained <- glm(yy ~ . , family="binomial",data=mydata, control = list(maxit = 1000))			
	
	trained$feat=apTrain$pathways.selected
	trained$crog=apTrain$selectedGenes	
	
	trained$graph = graph.adjacency(Gsub[trained$crog,trained$crog], mode="undirected")
	
	res =list(trained=trained,apTrain=apTrain)
	return (res)
}


predictPac <- function(train=train,x=x, y=y,int=int, DEBUG=FALSE)
{	
	apTrain<-train$apTrain
	apTest=probeset2pathwayTst(x=x,apTrain=apTrain)	
	
	mydata1 = data.frame(ap=apTest)
	trained<-train$trained
	test    <-  predict(trained, newdata=mydata1, type="response") #response , #cat("\t\t test done!\n")
	
	label <- sign(as.numeric(y) - 1.5) # because y is a factor
	nn=length(y)	
	pred=test-1.5
	auc   <- calc.auc(pred,label)
		
	res=list(auc=auc,pred=test)	
	return(res)
}

######
## sum the graph set with same pathway
# for selection the top genes in training data
probeset2pathwayTrain <-function(x=x, y=y, int=int)
{
	#require(AnnotationDbi)
	#::
	x<-scale(x)
		
	# for (reactome.db)	
	#KEGGids = mget(int,reactomeKEGGEXTID2PATHID,ifnotfound=NA)
	
	#for KEGG.db
	KEGGids = mget(int,KEGGEXTID2PATHID,ifnotfound=NA)
	KEGGids = KEGGids[which(KEGGids!="NA")]
	##
	kidList <- unlist(KEGGids)
	kidList <- kidList[which(kidList!="NA")]		
	ukidList <-unique(kidList)
	
	## form probset belond to KeegID, 	
	nCol = length(ukidList)	
	nRow = nrow(x)
	nColx = ncol(x)
	ap <- matrix(0, ncol=nCol, nrow=nRow)
	
	
	i=1
	selected.genes = list()
	netscore.list = list()
	pathway.list = list()
	for(p in ukidList)	
	{	#cat("\n at i:",i,"\n")
		index = unlist(which(kidList==p))
		index=as.character(intersect(colnames(x),index))
		len.k=length(index)	
		kk=1
		if(len.k==0)
		{
			selected.genes[[i]] = c()
			ap[,i] = 0
		}
		else if(len.k >= 1)
		{			
			tTest= matrix(0,len.k)
			
			for(ii in 1: len.k)
				tTest[ii] = abs(t.test(x=x[y==1, index[ii]], y=x[y==-1, index[ii]])$statistic)				

			tTest.order = order(tTest, decreasing=T)		
			netscore.new=0.0
			if(len.k>=2)
			{		
				tTest.order = order(tTest, decreasing=T)				
				app.tmp = x[, index[tTest.order[1]]]
				netscore = abs(t.test(x=app.tmp[y==1], y=app.tmp[y==-1])$statistic)	
				for(kk in 2: len.k)			
				{
					indG  = index[tTest.order[1:kk]]	
					app.tmp = rowSums(x[, indG])/sqrt(length(indG))			
					netscore.new = abs(t.test(x=app.tmp[y==1], y=app.tmp[y==-1])$statistic)						
					if(netscore.new > netscore)
					{	netscore = netscore.new		}
					else
					{	break}
				}
							
			}				
			indC=index[tTest.order[1:kk]]
			selected.genes[[i]] =  index[tTest.order[1:kk]]
			netscore.list[[i]] = netscore.new
			pathway.list[[i]] = p
			if(length(indC)>1)
			{ap[,i] = rowSums(x[, indC])/sqrt(length(indC))				}
			else 
			{ap[,i] = x[, indC]}
			#cat("\n ap:", head(ap[,i]), "\tat i:",i,
		}
		#cat(" #(selection genes)=",length(index[tTest.order[1:kk]]),"\nhead (ap):",head(ap[,i]),"at i:",i,"\n")	
		i=i+1
	}
	
	ap2=ap[,which(colSums(ap)!=0)]	
	
	netscoreMatrix=unlist(netscore.list)
	top.id <- which(netscoreMatrix>quantile(netscoreMatrix, .5))
	pathway.top = unique(unlist(pathway.list)[top.id])
	selected.gene.top=unlist(selected.genes[top.id])
		
	ap.top=ap2[, top.id]
	colnames(ap.top)=pathway.top 
	rownames(ap)= rownames(x)
	topList=list(path=pathway.top, crogs=selected.gene.top)	
	#cat("path ids:", pathway.top,"\n")
	res = list(ap= ap.top, selectedGenes = selected.gene.top, scaling=attributes(x), pathways.selected=pathway.top, top.list=topList)
	return (res)
}

#for test data
probeset2pathwayTst <-function(x=x,apTrain=apTrain)
{	
	# take selected pathways
	# take selected genes within pathways
	# apply z-transformation (with parameters adapted during training phase)
	# compute pathway activities
	x=scale(x, center=apTrain$scaling$"scaled:center", scale=apTrain$scaling$"scaled:scale")
	#cat("sacle\n")
	pathways=apTrain$top.list$path
	crogs=apTrain$top.list$crogs
	#kids = mget(selectedGenes,KEGGEXTID2PATHID, ifnotfound=NA)
	pid=  unlist(pathways)
	#cat("path ids:", pid,"\n")
	#Gids = mget(pid,KEGGPATHID2EXTID,ifnotfound=NA)		
	## form probset belond to KeegID, 	
	nCol = length(pid)
	nRow = nrow(x)
	nColx = ncol(x)	
	ap <- matrix(0, ncol=nCol, nrow=nRow)			
	
	for(i in 1:nCol)	
	{	#cat("\n at i:",i,"\n")
		index = crogs[[i]]
		index=as.character(index)
		if(length(index)== 1) ap[,i] = x[, index]
		else if(length(index)>= 2) ap[,i] = rowSums(x[, index])/sqrt(length(index))		
		#pid[i]=paste("ap.",pid[i],sep="")			
	}		
	#cat(" ap", "nCol:",nCol,"\tlegnth(pid):", length(pid),"\n")
	colnames(ap) = pid
	#cat("ap+pathway\n")
	rownames(ap) <- rownames(x)
	return (ap)
}
#apTst=probeset2pathwayTst(x=xtst,selectedGenes=apTrn$selectedGenes, pathways.selected=apTrn$pathways.selected,scaling.attrs=apTrn$scaling.attrs)

