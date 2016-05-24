evallist <- list()


evallist[[1]] <- function(tree,data,parms){
#mse for cart
    yname <- as.character(tree$terms[[2]])
	pred <- predict(tree,data)
	actual <- data[,yname]
	sse <- sum( (pred-actual)^2)
	#sst <- sum((actual - mean(actual))^2)
	#r2 <- 1 - (sse/sst)
	#return( r2 )
	return( sse/length(pred) )
}

evallist[[2]] <- function(tree,data,parms){
#miscl rate for cart
    yname <- as.character(tree$terms[[2]])
	pred <- predict(tree,data,type="class")
	actual <- data[,yname]
	return( sum(pred!=actual)/length(pred))
}


evallist[[3]]  <- function(tree,newdata,parms){
# purity for classification, evaluate by gini
    yname <- as.character(tree$terms[[2]])
	pred <- predict(tree,data,type="class")
	actual <- data[,yname]
}

names(evallist) <- c("anova","class","class_purity")


xval.OneSidedvsCART <- function(tree.obj,data,numxval,insamp.frac=.5,eval.fcn.list,nsfrac=TRUE)
{
#for doing xvals of one-sided methods versus non one-sided.
# if nsfrac=TRUE, then we evaluate the nodesize on tree.obj in context
# of original data and fix the fraction for the xvals. If false, we use the
# absolute number.
#
# For high means/low means + classification we need a class of interest. Arbitrarily,
# we choose the first 3 listed in levels(y) as a sample. 
#
	if(!inherits(tree.obj, "itree")) stop("Not a legitimate itree object")
	treemethod <- tree.obj$method
	class_of_interest <- NA
	mm <- tree.obj$method
	if(mm %in% c("class","class_extremes","class_purity")){
		mm <- "class"
		ll <- unique(data[,1]) #counting on the fact that y is in first column
		if(length(ll)>2){class_of_interest <- as.character(ll[1:3])}
		else{
			class_of_interest <- c(as.character(ll),NA)
		}
	}

	if(mm %in% c("anova_extremes","regression_purity","regression_extremes")){
		mm <- "anova"
	}
	n <- length(tree.obj$y)
	p <- dim(data)[2]
	subsize <- round(n*insamp.frac)
	idx <- 1:n
	ns = eval(tree.obj$call$minbucket)
	if(nsfrac==TRUE){ #reset to proper fraction
			ns <-  round(ns*insamp.frac)	
	}

	#do a match to get correct evaluation function for a given method
	evalmethod <- eval.fcn.list[[match(mm, names(eval.fcn.list))]] 

	stats <- matrix(0,ncol=5,nrow=numxval)
	for(i in 1:numxval){
		idx.sub1 <- sample(idx,size=subsize,replace=TRUE)
		idx.sub2 <- setdiff(idx,idx.sub1) #those not in 1st sumsample
		insample <- data[idx.sub1,]; outsample <- data[idx.sub2,]
		
		#do regular cart
		regularcart <- itree(eval(tree.obj$call$formula),data=insample,method=mm,
						minbucket= ns,
						minsplit= ns,
						cp = eval(tree.obj$call$cp),
						parms = tree.obj$parms,xval=0)
		
		stats[i,1] <- evalmethod(regularcart,outsample,tree.obj$call$parms)
		

	 	if(mm=="anova"){
	 	
	 		#do one-sided
			purity <- itree(eval(tree.obj$call$formula),data=insample,method="purity",
						minbucket= ns,
						minsplit= ns,
						cp = eval(tree.obj$call$cp),
						parms = tree.obj$parms,xval=0)
			stats[i,2] <- evalmethod(purity,outsample,tree.obj$call$parms)
	 	
	 		#do high means
			exthm <- itree(eval(tree.obj$call$formula),data=insample,method="extremes",
							minbucket= ns,
							minsplit= ns,
							cp = eval(tree.obj$call$cp),
							parms = 1,xval=0)	
			stats[i,3] <- evalmethod(exthm,outsample,tree.obj$call$parms)
									
			extlm <- itree(eval(tree.obj$call$formula),data=insample,method="extremes",
							minbucket= ns,
							minsplit= ns,
							cp = eval(tree.obj$call$cp),
							parms = -1,xval=0)	
			stats[i,4] <- evalmethod(extlm,outsample,tree.obj$call$parms)
		}
		else{ #need to specify class of interest.
		
			 #do one-sided
			purity <- itree(eval(tree.obj$call$formula),data=insample,method="purity",
						minbucket= ns,
						minsplit= ns,
						cp = eval(tree.obj$call$cp),xval=0)
			stats[i,2] <- evalmethod(purity,outsample,tree.obj$call$parms)

			interest1 <- itree(eval(tree.obj$call$formula),data=insample,method="extremes",
							minbucket= ns,
							minsplit= ns,
							cp = eval(tree.obj$call$cp),
							parms = list(classOfInterest=class_of_interest[1]),xval=0)	
			stats[i,3] <- evalmethod(interest1,outsample,tree.obj$call$parms)
			 
			interest2 <- itree(eval(tree.obj$call$formula),data=insample,method="extremes",
							minbucket= ns,
							minsplit= ns,
							cp = eval(tree.obj$call$cp),
							parms = list(classOfInterest=class_of_interest[2]),xval=0)
			stats[i,4] <- evalmethod(interest2,outsample,tree.obj$call$parms)				
			if(!is.na(class_of_interest[3])){

				interest3 <- itree(eval(tree.obj$call$formula),data=insample,method="extremes",
							minbucket= ns,
							minsplit= ns,
							cp = eval(tree.obj$call$cp),
							parms = list(classOfInterest=class_of_interest[3]),xval=0)
				stats[i,5] <- evalmethod(interest3,outsample,tree.obj$call$parms)
			} 
		
		}
	}
	#return
	if(mm=="anova"){
		stats <- stats[,1:4]
		colnames(stats) <- c("cart","purity","high.means","low.means")
	}
	if(mm != "anova"){
		if(is.na(class_of_interest[3])){
			colnames(stats) <- c("cart","purity",class_of_interest[1:2],"")
		}
		else{
			colnames(stats) <- c("cart","purity",class_of_interest)
		}
	}
	return(as.data.frame(stats))
}

xvalpenalty <- function(tree.obj,data,penaltytype,numxval,penalty.vec,
				eval.fcn.list){
	#assume target is called y
	# assumes interp_param1 is not specified
	#constants
	mm <- tree.obj$method
	if(mm %in% c("class_extremes","class_purity")){
		mm <- "class"
	}
	if(mm %in% c("anova_extremes","regression_purity","regression_extremes")){
		mm <- "anova"
	}
	n <- length(tree.obj$y)
	p <- dim(data)[2]
	subsize <- round(n * (2/3))
	idx <- 1:n
		
	#do a match to get correct evaluation function
	#for a given method
	evalmethod <- eval.fcn.list[[match(mm, names(eval.fcn.list))]] 

	stats <- matrix(0,ncol=length(penalty.vec),nrow=numxval)
	unpenalized <- rep(0,numxval)
	for(i in 1:numxval){
		idx.sub1 <- sample(idx,size=subsize,replace=TRUE)
		idx.sub2 <- setdiff(idx,idx.sub1) #those not in 1st sumsample
		insample <- data[idx.sub1,]; outsample <- data[idx.sub2,]
		
		#first do unpenalized
		temp <- itree(tree.obj$call$formula,data=insample,method=mm,
						minbucket= eval(tree.obj$call$minbucket),
						minsplit= eval(tree.obj$call$minsplit),
						cp = eval(tree.obj$call$cp),
						parms = tree.obj$call$parms,xval=0)
		xx <- evalmethod(temp,outsample,tree.obj$call$parms)
		unpenalized[i] <- xx
			
		#iterate thru penalty vals
		for(j in 1:length(penalty.vec)){
		
			#fit a tree with penalty j to i-th subset
			temp <- itree(tree.obj$call$formula,data=insample,method=mm,
						minbucket= eval(tree.obj$call$minbucket),
						minsplit= eval(tree.obj$call$minsplit),
						cp = eval(tree.obj$call$cp),
						parms = tree.obj$call$parms,
						xval=0,
						interp_param1 = penalty.vec[j],penalty=penaltytype)
			xx <- evalmethod(temp,outsample,tree.obj$call$parms) #function of fitted model and out-of-sample data
			stats[i,j] <- xx
			
		}#end loop thru penalties
	}#end loop thru bootstrap samples
	
	#clean up 
	cnames <- paste("k=",penalty.vec,sep="")
	rnames <- paste("sample",1:numxval,sep="")
	colnames(stats) <- cnames; rownames(stats) <- rnames
	names(unpenalized) <- rnames
	return(list(unpenalized=unpenalized,xvals=stats,penalty=penalty.vec))
}

xval.withWithoutPenalty <- function(tree.obj,data,penaltytypes,penalty.vec,numxval,insamp.frac,eval.fcn.list,nsfrac=TRUE){
	
	penalty.vec <- sort(penalty.vec) #increasing
	treemethod <- tree.obj$method
	class_of_interest <- NULL
	mm <- tree.obj$method
	if(mm %in% c("class","class_extremes","class_purity")){
		mm <- "class"
	}
	if(mm %in% c("anova_extremes","regression_purity","regression_extremes")){
		mm <- "anova"
	}

	n <- length(tree.obj$y)
	p <- dim(data)[2]-1
	subsize <- round(n*insamp.frac)
	idx <- 1:n
	ns = eval(tree.obj$call$minbucket)
	if(nsfrac==TRUE){ #reset to proper fraction
		ns <-  round(ns*insamp.frac)	
	}

	#do a match to get correct evaluation function for a given method
	evalmethod <- eval.fcn.list[[match(mm, names(eval.fcn.list))]] 
	
	holdout.stats <- matrix(0,nrow=numxval,ncol=(length(penaltytypes)+1))
	insample.stats <- matrix(0,nrow=numxval,ncol=(length(penaltytypes)+1))
	bestk.vals <- matrix(0,nrow=numxval,ncol=length(penaltytypes))
	
	########### BEGIN LOOP THRU XVALS
	for(i in 1:numxval){
		# BOOTSTRAP
		idx.sub1 <- sample(idx,size=subsize,replace=TRUE)
		idx.sub2 <- setdiff(idx,idx.sub1) #those not in 1st sumsample
		insample <- data[idx.sub1,]; outsample <- data[idx.sub2,]

		#fit the nonpenalized tree.
		if(mm=="anova"){
			ppp <- tree.obj$parms
		}
		else{
			ppp <- tree.obj$call$parms
			if(treemethod=="class_extremes"){
				ppp <- eval(ppp)
			}
		}
		nonpenalized <- itree(eval(tree.obj$call$formula),data=insample,method=treemethod,
					minbucket=ns,
					minsplit= ns,
					cp = eval(tree.obj$call$cp),
					parms = ppp,xval=0)
		
		holdout.stats[i,1]  <- evalmethod(nonpenalized,outsample,NULL)
		insample.stats[i,1] <- evalmethod(nonpenalized,insample,NULL)

		performance.vs.penalty1 <- matrix(0,nrow=length(penalty.vec),ncol=length(penaltytypes))
		performance.vs.penalty2 <- matrix(0,nrow=length(penalty.vec),ncol=length(penaltytypes))
		#ITERATE THRU EACH PENALTY VALUE
		for(k in 1:length(penalty.vec)){
			#ITERATE THRU EACH PENALTY METHOD, evaluating in and out of sample
			for(ptype in 1:length(penaltytypes)){
				temp <- itree(eval(tree.obj$call$formula),data=insample,method=treemethod,
						minbucket= ns,
						minsplit= ns,
						cp = eval(tree.obj$call$cp),
						parms = ppp,
						xval=0,
						interp_param1 = penalty.vec[k],penalty=penaltytypes[ptype])
				#plot(temp); text(temp)
				performance.vs.penalty1[k,ptype] <- evalmethod(temp,insample,NULL)
				performance.vs.penalty2[k,ptype] <- evalmethod(temp,outsample,NULL)
			}#END LOOP THRU EACH PENALTY METHOD			
		}#END LOOP THRU PENALTY VALS
		
		#FOR EACH PENALTY METHOD CHOOSE THE BEST K
		target <- insample.stats[i,1] #performance of nonpenalized
		print(target)
		acceptable <- (performance.vs.penalty1 <= (1.1*target))
		
		#choose highest penalty with acceptable performance 
		#recall penalty.vec is  incr
		bestk <- apply(apply(acceptable,2,cumsum)*acceptable,2,which.max)
		#RECORD SOME STUFF
		#print(performance.vs.penalty1)
		bestk.vals[i,] <- penalty.vec[bestk]
		for(p in 1:length(penaltytypes)){
			holdout.stats[i,1+p] <- performance.vs.penalty2[bestk[p],p]				
			insample.stats[i,1+p] <- performance.vs.penalty1[bestk[p],p]				
		}
		print(paste("xval ",i,":"," bestk=",bestk[1]," ",bestk[2],sep=""))
	}#MOVE ON TO THE NEXT XVAL

	#clean up the data and return it.	
	colnames(bestk.vals) <- penaltytypes
	rownames(bestk.vals) <- rownames(insample.stats) <- rownames(holdout.stats) <- paste("xval",(1:numxval),sep="")
	colnames(insample.stats) <- colnames(holdout.stats) <- c(treemethod,paste(treemethod,penaltytypes,sep="."))
	list(insample=insample.stats,holdout=holdout.stats,kstar=bestk.vals,penalty.vec=penalty.vec)
}


plotxvals <- function(xvalobject){
	nxval <- length(xvalobject$unpenalized)
	pp <- xvalobject$penalty
	np <- mean(xvalobject$unpenalized)
	np.se <- sd(xvalobject$unpenalized)/sqrt(nxval)
	penalized.mean <- as.vector(apply(xvalobject$xvals,2,mean))
	penalized.se <- as.vector(apply(xvalobject$xvals,2,sd))/sqrt(nxval)
	ymax <- 2*max(penalized.se)+max(penalized.mean)
	ymax <- max(ymax,np+2*np.se)
	ymin <- min(penalized.mean) - 2*max(penalized.se)
	ymin <- min(ymin, np-2*np.se) 
	plot(x=pp,y=penalized.mean, ylim=c(ymin,ymax),type="b",
		main="CV stat vs penalty",xlab="penalty",ylab="cv",lwd=1.5)
	par(new=TRUE)
	plot(x=pp,y=penalized.mean+2*penalized.se, ylim=c(ymin,ymax),type="b",col="RED",axes=F,
		xlab='',ylab='',lty=2)
	par(new=TRUE)
	plot(x=pp,y=penalized.mean-2*penalized.se, ylim=c(ymin,ymax),type="b",col="RED",axes=F,
		xlab='',ylab='',lty=2)
	abline(b=0,a=np,lwd=2); abline(b=0,a=np+2*np.se,lty=2); abline(b=0,a=np-2*np.se,lty=2) 
}
