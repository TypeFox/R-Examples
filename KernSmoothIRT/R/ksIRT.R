ksIRT <-
function(responses, key, format, kernel = c("gaussian","quadratic","uniform"), itemlabels, weights, 
miss = c("option","omit","random.multinom","random.unif"), NAweight = 0, evalpoints, nevalpoints, 
bandwidth = c("Silverman","CV"), RankFun = "sum", SubRank, thetadist = list("norm",0,1), groups = FALSE){

	if(missing(key)) key <- NULL
	if(missing(itemlabels)) itemlabels <- NULL
	if(missing(weights)) weights <- NULL
	if(missing(evalpoints)) evalpoints <- NULL
	if(missing(SubRank)) SubRank <- NULL
	if(missing(nevalpoints)) nevalpoints <- NULL
	if(missing(format)) format <- NULL
	

	
	if(is.null(evalpoints) & is.null(nevalpoints)){nevalpoints <- 51}
	else if(is.null(nevalpoints)){nevalpoints <- length(evalpoints)}
	
	if(class(responses) == "list") responses = do.call(cbind,responses)
	responses <- as.matrix(responses)
	
	kernel <- match.arg(arg = kernel, choices = c("gaussian","quadratic","uniform"))
	miss   <- match.arg(arg = miss, choices = c("option","omit","random.multinom","random.unif"))
	bandwidth <- match.arg(arg = bandwidth, choices = c("Silverman","CV"))
	RankFun <- match.fun(RankFun)


	if(miss=="omit"){
	## Remove subjects with missing responses.
        if(groups[1] != FALSE){ groups <- groups[complete.cases(responses)]}
		responses <- na.omit(responses)
	}

	nsubj <- nrow(responses);
	nitem<-ncol(responses);
	
	if(length(key) == 1){key <- rep(key, nitem)}

	Check<-Inputcheck(responses,key,format,itemlabels,weights,miss,evalpoints,bandwidth,nitem,nsubj,kernel,NAweight,nevalpoints,thetadist,groups,SubRank)
	## Make sure input variables are specified correctly or return a message
	if(Check!="NoErr"){stop(Check)}


	## Determine type of response for each item.
	if(length(format)==1){
		if(format[1]==1) scale<-rep(1,nitem)
		else if(format[1]==2)scale<-rep(0,nitem)
		else if(format[1]==3){
			scale<-rep(3,nitem) 
			print("Only OCC plots are available for nominal items. Use option axis = 'distribution' when plotting.")
		}
	}
	else{scale <- format; scale[scale == 2] <- 0;}
	
	if(is.null(itemlabels)){itemlabels<-colnames(responses)}
	## Default item labels

	responses<-apply(responses,2,handlemiss,miss=miss)
	## This function (handlemiss) takes a single subjects' responses and treats missing values accordingly.
	optsbyitem<-list()


	#########################################################################################################################
	## The next several lines take the response matrix, and turn it into a larger binary matrix with three additional columns. 
	## The first, is the item number, the second is the option number, the third is the option weight.
	## The result is a matrix with the number of columns equal to the number of subjects plus 3.
	## The number of rows is equal to the total number of options on the exam.
	## Other than the first 3 columns, each entry in the matrix is either a 0 or 1, indicating whether or not the subject selected
	## that particular option on that particular item.
	##########################################################################################################################

	for(i in 1:ncol(responses)){
		optsbyitem[[i]] <- unique(responses[,i])
	}

	fullresponses<-matrix(0,length(unlist(optsbyitem)),ncol=3+nsubj)

		crow<-0
		for(i in 1:nitem){
			for(j in 1:length(optsbyitem[[i]])){
				crow<-crow+1

					fullresponses[crow,1:2]<-c(i,optsbyitem[[c(i,j)]])
					
					if(!is.null(weights)){	
						fullresponses[crow,3]<-getweight(item=i,option=fullresponses[crow,2],weights=weights[[i]],NAweight=NAweight)
					}
					if(!is.null(key)){
			
						fullresponses[crow,3]<-getweight(item=i,option=fullresponses[crow,2],scale=scale[i],key=key[i],NAweight=NAweight)}						

				
				}


		}


	fullresponses0<-fullresponses
	optitwgtresp<-make_mat(A=fullresponses,B=responses)

	################################################################################
	################################################################################


	 APPLYFUN <- function(col,function_name,weights){
		FUN <- match.fun(function_name)
		FUN(weights * col)
	}
	 
	ncorrectex<-apply(optitwgtresp[,-c(1:3)],2,APPLYFUN,"sum",optitwgtresp[,3])

	## Turn the first argument of thetadist (the distribution ex. norm, unif, beta) into the quantile function.
	qdist <- get(paste("q",thetadist[[1]],sep=""),mode="function")
	## Rank order each eqaminee according to the number correct.
	if(is.null(SubRank)) {
		rankscores<-rank(apply(optitwgtresp[,-c(1:3)],2,APPLYFUN,RankFun,optitwgtresp[,3]),ties.method="first")
	} 
	else{
		rankscores <- SubRank
	}
	## Turn that rank into a quantile of the selected proability distribution with parameters as specified in the other arguments of the thetadist option.
	subjtheta<-eval(parse(text=(paste("qdist((rankscores)/(nsubj+1),",paste(unlist(thetadist[-1]),collapse=","),")",sep=""))))

	## For plotting purposes, find the appropriate subjscoresummary.
	quantsevalpoints<-quantile(subjtheta,probs=c(.05,.25,.50,.75,.95))

	

	if(is.null(evalpoints)){
			## evalpoints is a vector at which we evaluate the kernel smoother and plot 
			## if not specified, then use equal spacing.
			lim1<-qdist(1/(nsubj+1),thetadist[[2]],thetadist[[3]])
			lim2<-qdist(nsubj/(nsubj+1),thetadist[[2]],thetadist[[3]])
			evalpoints<-seq(from=lim1, to=lim2, length.out=nevalpoints);

		}
		

	
	## Find the closest actual scores to each of the evalpoints points.
	scoresatevalpoints<-quantile(ncorrectex,(1:nevalpoints)/nevalpoints)
	


	if(bandwidth[1]=="CV"){		
		## If we are using cross-validation to find the optimal bandwidth for each item
		h<-numeric()
		h0<-numeric()
		torep<-table(optitwgtresp[,1])

		for(i in 1:nitem){
			
			mat<-optitwgtresp[which(optitwgtresp[,1]==i),-c(1:3)]	
			print(paste("Determining Smoothing Bandwidth Item: ",i,sep=""))
			## This is computationally intensive task, so output a message for the impatient.
			h0[i]<-CrossV(answered0=mat,subjtheta0=subjtheta,kernel0=kernel)
			## This function returns the optimal bandwidth for each item. 
			h<-c(h,rep(h0[i],torep[i]))

		}
	}

	else if(bandwidth[1]=="Silverman"){
		## If we are not using cross-validation and haven't specified a bandwidth, then we 
		## use the default option which is h (below) for all items.
		if(thetadist[[1]]=="norm"){
			sighat <- thetadist[[3]]			
			
		}
		else{
			sighat <- sd(subjtheta)
		}		
			h<-rep(1.06*sighat*nsubj^(-.2),nrow(optitwgtresp))
			h0<-rep(h[1],nitem)
	}

	else{
		## If we have specified a bandwidth for each item, then use it on each option 
		## within each item.
		torep<-table(optitwgtresp[,1])
		h0<-bandwidth;
		h<-numeric()

		for(i in 1:nitem){
			h<-c(h,rep(h0[i],torep[i]))
		}

	}



	## Here perform the actual smoothing. 
	ICC<-matrix(0,nrow=nrow(optitwgtresp),ncol=nevalpoints)
	kernelweights<-matrix(0,nrow=nrow(optitwgtresp),ncol=nevalpoints)
	stderr<-matrix(0,nrow=nrow(optitwgtresp),ncol=nevalpoints)

	if(kernel=="gaussian"){ktog<-1;}
	if(kernel=="quadratic"){ktog<-2;}
	if(kernel=="uniform"){ktog<-3;}


	for (i in 1:nrow(optitwgtresp)){
		## C++ function that does the smoothing. 
		retval<-smoother(A = h[i], B = subjtheta, C = evalpoints, D = optitwgtresp[i,-c(1:3)], E = ktog) 
		## Return the standard error matrix, the option characteristic matrix (inappropriately named ICC) and the smoothing weight matrix
		ICC[i,]<-retval[["ICC"]]
		stderr[i,]<-retval[["stderr"]]
		kernelweights[i,]<-retval[["weights"]]
	}

	
	

	## Attach item, option number and weight to smoothed OCC matrix and standard error matrix.
	Probs<-cbind(optitwgtresp[,c(1:3)],ICC)
	stderrs<-cbind(optitwgtresp[,c(1:3)],stderr)
	stderrs[which(is.na(stderrs))]<-0

	## Perform DIF analysis if specified.
	if(groups[1]==FALSE){DIF<-FALSE; grps<-NULL;}
	else{
		grps<-unique(groups)
		## Recursive call for each of the groups.
		DIF<-lapply(grps,function(x)ksIRT(responses=responses[which(groups==x),],key=key,format=format,kernel=kernel,itemlabels=itemlabels,weights=weights,miss=miss,NAweight=NAweight,evalpoints=evalpoints,nevalpoints=nevalpoints,bandwidth=bandwidth,thetadist=thetadist,groups=FALSE))

	}
	
	
	evals <-  Probs[,3] %*% Probs[,-c(1:3)]
	quantsex<-quantile(ncorrectex,probs=c(.05,.25,.50,.75,.95))

	
	
	## Get each item score
	this<-cbind(optitwgtresp[,1],t(apply(optitwgtresp,1,function(xxx)xxx[-c(1:3)]*xxx[3])))
	corr<-numeric(nitem)


	
	Probs0 <- numeric(0)
	for(i in 1:nitem){
		thisit<-which(this[,1]==i)
	
				summed<-apply(this[thisit,-1],2,sum)	
				if(scale[i] != 3) corr[i] <- cor(ncorrectex,summed,use="complete.obs")	
				else corr[i] <- 0
				Probs00 <- cbind(Probs[Probs[,1]==i,1:3],apply(Probs[Probs[,1]==i,-c(1:3)],2,function(x)x/sum(x)))				
				
		
		## Find the biserial correlation between items and overall score.
		
	
	
		
		Probs0<-rbind(Probs0,Probs00)
	
	
	
	}


	LIK <- list(nsubj)	
	MLE <- numeric(nsubj)
	MLEthet <- numeric(nsubj)
	
	for(i in 1:nsubj){
		## Get the likelihood and MLE for each subject.
		responses <- optitwgtresp[, i + 3]
		GOOD <- Probs0[,-c(1:3)] * responses
		GOOD[GOOD==0] <- 1
		
		LIK[[i]] <- apply(GOOD,2,prod)
		
		LIK[[i]] <- LIK[[i]]/max(LIK[[i]])
		MLE[i] <- evals[which.max(LIK[[i]])]
		MLEthet[i] <- evalpoints[which.max(LIK[[i]])]
	}

	## return values.
	toret<-list(binaryresp=optitwgtresp,OCC=Probs0,stderrs=stderrs,subjscore=ncorrectex,itemlabels=itemlabels,evalpoints=evalpoints, subjscoresummary=quantsex,subjthetasummary=quantsevalpoints,kernelweights=kernelweights,scale=scale,thetadist=thetadist, subjtheta=subjtheta,bandwidth=h0,nitem=nitem,nsubj=nsubj,nevalpoints=nevalpoints,DIF=DIF,groups=grps,expectedscores=evals,itemcor=corr, subjscoreML = MLE, subjthetaML = MLEthet, RCC = LIK, format=format)


	class(toret)<-"ksIRT"

	return(toret)

}
