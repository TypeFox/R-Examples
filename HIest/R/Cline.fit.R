Cline.fit <- 
function(Data,By=NULL,S=NULL,model,Start=NULL,Methods=NULL,iterations=99,SD=NULL,headstart=TRUE,Grid=FALSE,ploidy=2,trim=0,include=1:ncol(Data)){
# ploidy can (not) be a vector specifying different ploidy for each locus
	models <- rep(FALSE,6)
	Models <- c("multinom","binom","logit.logistic","Barton","Beta","Richards")
	if(sum(is.na(match(model,Models)))>0) stop("models must be one or more of: \n logit.logistic, Barton, Beta, Richards, multinom, binom\n")
	
	n.mods <- length(model)
	n <- dim(Data)[1]
	L <- dim(Data)[2]
	
	if(!is.null(By)) { #aggregate data by the specified factor
		mfun <- function(x,ploidy,na.rm=TRUE,trim=0,include){mean((x/ploidy)[include],na.rm=TRUE,trim=trim)}
		nfun <- function(x,ploidy){sum(!is.na(x))*ploidy}
		X.S <- apply(aggregate(Data,by=list(By),FUN=mfun,ploidy=ploidy,na.rm=TRUE,trim=trim)[-1],1,mean,na.rm=TRUE)
		X.n <- aggregate(Data,by=list(By),FUN=nfun,ploidy=ploidy)[,-1]
		X.x <- aggregate(Data,by=list(By),FUN=sum,na.rm=TRUE)[,-1]
		X.h <- aggregate(Data==1,by=list(By),FUN=sum,na.rm=TRUE)[,-1]
		X.0 <- aggregate(Data==0,by=list(By),FUN=sum,na.rm=TRUE)[,-1]
		X.2 <- aggregate(Data==2,by=list(By),FUN=sum,na.rm=TRUE)[,-1]
}
	if(is.null(By)){
		mfun <- function(x,ploidy,na.rm=TRUE,trim=0,include){
			mean((x/ploidy)[include],na.rm=na.rm)
			}
		X.S <- apply(Data,1,mfun,ploidy=ploidy,na.rm=TRUE,include=include)
		X.n <- matrix(ploidy,nrow=n,ncol=L,byrow=TRUE)
		X.x <- Data
		X.h <- Data==1
		}
	if(!is.null(S)){X.S <- S}
	if(length(X.S)!=dim(X.x)[1]){return(expression("bad S"))}
	if(sum(is.na(X.S)>0)){
		cat("\nWarning: some observations with all NA entries are being omitted\n")
		X.out <- !is.na(X.S)
		X.S <- X.S[X.out]
		X.n <- X.n[X.out,]
		X.x <- X.x[X.out,]
		}
		
	if(is.null(Start)){
		p1 <- 1e+5
		Start=list(llogit=c(u=0,v=1),bart=c(a=0,b=0),beta=c(.5,2),rich=c(p1,1-p1,-2*log((p1-1)/p1),1/2))
		}
	if(is.null(Methods)){ Methods=list(logit.logistic="L-BFGS-B",Barton="L-BFGS-B",Beta="L-BFGS-B",Richards="L-BFGS-B") } # if Methods is not null, it must be a named list
	if(is.null(SD)){ SD <- list(llogit = c(.1,.1),bart=c(.1,.1),beta=c(.1,.1),rich=rep(.1,4)) }
	if(is.matrix(Data)){Locus.names <- colnames(Data)}
	if(is.data.frame(Data)){Locus.names <- names(Data)}
	
	if(match(Models,model,nomatch=0)[1]>0){ # multinomial
		require(nnet)
		MN.out <- data.frame(int.1=NA,slp.1=NA,int.2=NA,slp.2=NA,lnL=NA,AICc=NA)
		models[1] <- TRUE}
	if(match(Models,model,nomatch=0)[2]>0){ # binomial
		require(stats)
		BN.out <- data.frame(a=NA,b=NA,lnL=NA,AICc=NA)
		models[2] <- TRUE}
	if(match(Models,model,nomatch=0)[3]>0){ # logit-logistic
		llogit.out <- data.frame(u=NA,v=NA,lnL=NA,AICc=NA)
		models[3] <- TRUE}	
	if(match(Models,model,nomatch=0)[4]>0){ # Barton
		bart.out <- data.frame(a=NA,b=NA,lnL=NA,AICc=NA)
		models[4] <- TRUE}	
	if(match(Models,model,nomatch=0)[5]>0){ # Beta
		beta.out <- data.frame(mu=NA,nu=NA,lnL=NA,AICc=NA)
		models[5] <- TRUE}	
	if(match(Models,model,nomatch=0)[6]>0){ # Richards
		rich.out <- data.frame(U=NA,L=NA,slope=NA,m=NA,lnL=NA,AICc=NA)
		models[6] <- TRUE}	
	best.fit <- character(); bf <- numeric()
	
	for(i in 1:L){
	if(length(unique(na.omit(Data[,i])))>1){
		if(models[1]){
			if(is.null(By)){y <- cbind(X.x[,i]==0,X.x[,i]==1,X.x[,i]==2)}
			if(!is.null(By)){y <- cbind(X.0[,i],X.h[,i],X.2[,i])}
			here <- X.n[,i]>0
		FIT <- gcline.fn(x=X.S[here],n=X.n[here,i],y=y[here,],start=NULL,model="multinom",method=NULL,iterations=iterations,SD=NULL,headstart=headstart)
			if(!is.null(dim(FIT$est))){ # if there are three genotypes
				MN.out[i,1:2] <- FIT$est[1,]
				MN.out[i,3:4] <- FIT$est[2,]
			}
			if(is.null(dim(FIT$est))){ # if there are only two genotypes
				cat("\nWarning -- Marker",colnames(Data)[i],"has only two genotypes for multinomial regression.\n")
				if(min(Data[,i],na.rm=TRUE)==0 & max(Data[,i],na.rm=TRUE)==1){
					MN.out[i,1:2] <- FIT$est
					MN.out[i,3:4] <- c(-25,0)
				}
				if(min(Data[,i],na.rm=TRUE)==1 & max(Data[,i],na.rm=TRUE)==2){
					MN.out[i,1:2] <- c(25,0)
					MN.out[i,3:4] <- FIT$est
				}
				if(min(Data[,i],na.rm=TRUE)==0 & max(Data[,i],na.rm=TRUE)==2){
					MN.out[i,1:2] <- c(-25,0)
					MN.out[i,3:4] <- FIT$est
				}
			}
		MN.out[i,5] <- FIT$lnL[1]; MN.out[i,6] <- bf[1] <- FIT$AICc
		}
		if(models[2]){
		FIT <- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=NULL,model="logistic",method=NULL,iterations=iterations,SD=NULL,headstart=headstart)
		BN.out[i,1:2] <- FIT$est; BN.out[i,3] <- FIT$lnL[1]; BN.out[i,4] <- bf[2] <- FIT$AICc
		}
		if(models[3]){
		FIT <- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=Start$llogit,model="logit-logit",method=Methods$logit.logistic,iterations=iterations,SD=SD$llogit,headstart=headstart)
		llogit.out[i,1:2] <- FIT$est; llogit.out[i,3] <- FIT$lnL[1]; llogit.out[i,4] <- bf[3] <- FIT$AICc
		}
		if(models[4]){
		FIT <- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=Start$bart,model="Barton",method=Methods$Barton,iterations=iterations,SD=SD$bart,headstart=headstart)
		bart.out[i,1:2] <- FIT$est; bart.out[i,3] <- FIT$lnL[1]; bart.out[i,4] <- bf[4] <- FIT$AICc
		}
		if(models[5]){
		FIT<- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=Start$beta,model="Beta",method=Methods$Beta,iterations=iterations,SD=SD$beta,headstart=headstart)
		beta.out[i,1:2] <- FIT$est; beta.out[i,3] <- FIT$lnL[1]; beta.out[i,4] <- bf[5] <- FIT$AICc
		}
		if(models[6]){
		FIT<- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=Start$rich,model="Richards",method=Methods$Richards,iterations=iterations,SD=SD$rich,headstart=headstart)
		rich.out[i,1:4] <- FIT$est; rich.out[i,5] <- FIT$lnL[1]; rich.out[i,6] <- bf[6] <- FIT$AICc
		}
		best.fit[i] <- Models[which.min(bf)]
	}
		if(length(unique(na.omit(Data[,i])))==1){
			cat("\nMarker",colnames(Data)[i],"is fixed in the sample.\n")
		if(unique(na.omit(Data[,i]))==0){
			if(models[1]) MN.out[i,] <- c(-25,0,-25,0,NA,NA)
			if(models[2]) BN.out[i,] <- c(-25,0,NA,NA)
			if(models[3]) llogit.out[i,] <- c(0,-25,NA,NA)
			if(models[4]) bart.out[i,] <- c(-25,-25,NA,NA)
			if(models[5]) beta.out[i,] <- c(1-.Machine$double.neg.eps,2,NA,NA)
			if(models[6]) rich.out[i,] <- c(1,1,0,0,NA,NA)
			best.fit[i] <- "fixed_0"
		}
		if(unique(na.omit(Data[,i]))==2){
			if(models[1]) MN.out[i,] <- c(-25,0,25,0,NA,NA)
			if(models[2]) BN.out[i,] <- c(25,0,NA,NA)
			if(models[3]) llogit.out[i,] <- c(0,25,NA,NA)
			if(models[4]) bart.out[i,] <- c(25,-25,NA,NA)
			if(models[5]) beta.out[i,] <- c(.Machine$double.xmin,2,NA,NA)
			if(models[6]) rich.out[i,] <- c(0,0,0,0,NA,NA)
			best.fit[i] <- "fixed_1"
		}
	}
	}
	outlist <- list()
	if(models[1]){rownames(MN.out)<-names(Data)
		D2 <- mahalanobis(MN.out[,1:4],colMeans(MN.out[,1:4]),cov(MN.out[,1:4]))
		P <- pchisq(D2,4,lower.tail=FALSE)
		outlier <- P < 0.05/L
		MN.out <- cbind(MN.out,D2=D2,P=P,outlier=outlier)
		outlist$multinom <- MN.out
	}
	if(models[2]){rownames(BN.out)<-names(Data)
		D2 <- mahalanobis(BN.out[,1:2],colMeans(BN.out[,1:2]),cov(BN.out[,1:2]))
		P <- pchisq(D2,2,lower.tail=FALSE)
		outlier <- P < 0.05/L
		BN.out <- cbind(BN.out,D2=D2,P=P,outlier=outlier)
		outlist$binom <- BN.out
	}
	if(models[3]){rownames(llogit.out)<-names(Data)
		D2 <- mahalanobis(llogit.out[,1:2],colMeans(llogit.out[,1:2]),cov(llogit.out[,1:2]))
		P <- pchisq(D2,2,lower.tail=FALSE)
		outlier <- P < 0.05/L
		llogit.out <- cbind(llogit.out,D2=D2,P=P,outlier=outlier)
		outlist$logit.logistic <- llogit.out
	}
	if(models[4]){rownames(bart.out)<-names(Data)
		D2 <- mahalanobis(bart.out[,1:2],colMeans(bart.out[,1:2]),cov(bart.out[,1:2]))
		P <- pchisq(D2,2,lower.tail=FALSE)
		outlier <- P < 0.05/L
		bart.out <- cbind(bart.out,D2=D2,P=P,outlier=outlier)
		outlist$Barton <- bart.out
	}
	if(models[5]){rownames(beta.out)<-names(Data)
		D2 <- mahalanobis(beta.out[,1:2],colMeans(beta.out[,1:2]),cov(beta.out[,1:2]))
		P <- pchisq(D2,2,lower.tail=FALSE)
		outlier <- P < 0.05/L
		beta.out <- cbind(beta.out,D2=D2,P=P,outlier=outlier)
		outlist$beta <- beta.out
	}
	if(models[6]){rownames(rich.out)<-names(Data)
		D2 <- mahalanobis(rich.out[,1:4],colMeans(rich.out[,1:4]),cov(rich.out[,1:4]))
		P <- pchisq(D2,4,lower.tail=FALSE)
		outlier <- P < 0.05/L
		rich.out <- cbind(rich.out,D2=D2,P=P,outlier=outlier)
		outlist$Richards <- rich.out
	}
	
	if(sum(models)>1){
		outlist$best.fit <- data.frame(marker=names(Data),best.fit=best.fit)
	}
	outlist
}