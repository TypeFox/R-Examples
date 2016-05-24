# Uses procedure outlined by Breiman to estimate the loss at each
# observation 1,...,N in 'data'.
# sampeFcn is a function that takes the inidices 1,..,N and returns
# a single sample to be the LEARNING data.  Default is size N bootstrap
# with replacement.
getOOBLoss <- function(model_tree.obj,data,nboot=100,
		sampleFcn = function(idx_vec){sample(idx_vec,replace=TRUE)},
		minsplit, minbucket,lossfcn)

{
	if(!inherits(model_tree.obj, "itree")) stop("Not legitimate itree object")
	
	#get the whole response vector.
	yname <- strsplit(deparse(model_tree.obj$call$formula),"~")[[1]][1]
	yname <- sub(" ","",yname)
	
	treemethod <- model_tree.obj$method  #passed to rpart
	mm <- model_tree.obj$method           #to figure out loss function
	if(mm %in% c("class_extremes","class_purity","class")){
		mm <- "class"
		Y <- match(data[,yname], levels(data[,yname])) #needs to be numeric, not string
		
	}
	if(mm %in% c("anova","regression_purity","regression_extremes")){
		mm <- "anova"
		Y <- data[,yname]
	}
	#print error for non anova/class methods
	if(!(mm %in% c("class","anova"))){
		stop("getOOBLoss not defined for this method.")
	}
	
	
	if(treemethod=="anova"){
		ppp <- model_tree.obj$parms
	}
	else{
		ppp <- model_tree.obj$call$parms
		if(treemethod=="class_extremes"){
			ppp <- eval(ppp)
		}
		if(treemethod=="regression_extremes"){
			ppp <- eval(ppp)
		}		
	}

	#some constants
	N <- nrow(data)
	p <- ncol(data)
	idx <- 1:N
	
	#what are the nodesizes?
	if(missing(minsplit)){ minsplit_final <- eval(model_tree.obj$control$minsplit) }
	else{
		if(is.numeric(minsplit)){
			if(minsplit<1){ minsplit_final <- round(N*minsplit)} #assume it's a fraction of N
			else{  minsplit_final <- round(minsplit)}
		}else{
			if(class(minsplit)!="function"){
				stop("Invalid minsplit argument. Pass a function of N, a fraction or an integer.")
			}
			minsplit_final <- minsplit(N)
		}
	}
	if(missing(minbucket)){minbucket_final <- eval(model_tree.obj$control$minbucket) }
	else{
		if(is.numeric(minbucket)){
			if(minbucket<1){ minbucket_final <- round(N*minbucket)} #assume it's a fraction of N
			else{  minbucket_final <- round(minbucket)}
		}else{
			if(class(minbucket)!="function"){
				stop("Invalid minbucket argument. Pass a function of N, a fraction or an integer.")
			} 
			minbucket_final <- minbucket(N) 
		}
	}
	
	# place to hold out-of-sample predictions 
	# holdout.predictions[i,j] = oob pred on obs of i in jth train sample
	# = NA if obs i is insample for run j.
	holdout.predictions <- matrix(NA,nrow=N,ncol=nboot)
	#bootstrap/xval runs...
	for(i in 1:nboot){
		idx.sub1 <- sampleFcn(idx) 
		idx.sub2 <- setdiff(idx,idx.sub1) #those not in 1st sumsample
		insample <- data[idx.sub1,]; outsample <- data[idx.sub2,]
		
		#get predictions
		temp <- itree(eval(model_tree.obj$call$formula),data=insample,method=treemethod,
				minbucket= minbucket_final,
				minsplit = minsplit_final,
				cp = eval(model_tree.obj$control$cp),
				parms = ppp,xval=0)
		#get predictions
		if(mm=="class"){
			preds <- predict(temp,outsample,type="class")
		}else{
			preds <- predict(temp,outsample)
		}
		holdout.predictions[idx.sub2,i] <- preds	
	}#end of bootstrap runs.
	
	#now clean up and format for output. 
	cnames <- paste("xval",(1:nboot),sep="")
	colnames(holdout.predictions) <- cnames
	
	#now assess mse, bias, variance
	num.not.na <- apply(holdout.predictions,1,function(temp){sum(!is.na(temp))})
	if(mm=="anova"){
		preds <- apply(holdout.predictions,1,function(temp){mean(temp[!is.na(temp)])})
		#varpred <- apply(holdout.predictions,1,function(temp){var(temp[!is.na(temp)])})
		YY <- matrix(rep(Y,nboot),nrow=N,ncol=nboot)
		if(missing(lossfcn)){
			YY <- (YY-holdout.predictions)^2
			YY[is.na(YY)] <- 0
			mses <- apply(YY,1,sum)/num.not.na
		}else{
			mses <- lossfcn(YY,preds)
		}
		return(list(bagpred=preds,holdout.predictions=holdout.predictions,avgOOBloss=mses))
	}
	else{ #classification
		Mode <- function(x) {
			x <- x[!is.na(x)]
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
		}
		preds <- apply(holdout.predictions,1,Mode)
		YY <- matrix(rep(Y,nboot),nrow=N,ncol=nboot)
		if(missing(lossfcn)){
			YY <- (YY != holdout.predictions)
			YY[is.na(YY)] <- 0
			miscl <- apply(YY,1,sum)/num.not.na
		}else{
			miscl <- lossfcn(YY,preds)
		}
		return(list(bagpred=preds,holdout.predictions=holdout.predictions,avgOOBloss=miscl))
	}
}		


estNodeRisk <- function(tree.obj,est_observation_loss){	
#row number in tree.obj$frame of the leaf node for each obs
	
	if(!inherits(tree.obj, "itree")) stop("Not legitimate itree object")
	
	ww <- tree.obj$where 
	nodes <- as.matrix(sort(unique(ww)),ncol=1)

	node.avg.loss <- function(nodenum){
		mean(est_observation_loss[ww==nodenum])
	} 
	node.sd.loss <- function(nodenum){
		sd(est_observation_loss[ww==nodenum])
	} 
		 
	avg.loss <- apply(nodes,1,node.avg.loss)
	sd.loss <- apply(nodes,1,node.sd.loss)

	temp <- list(est.risk=avg.loss,sd.loss=sd.loss)
	class(temp) <- "estNodeRisk"
	return(temp)	 
}