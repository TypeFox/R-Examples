

ROC.voting <- function(dat, grp, n = 200, healthy = NULL, seed = NULL, para1 = 200, para2 = 1, length = 40, display = FALSE, med = FALSE){

	para2 <- ceiling(para2)
	if(!is.null(seed)){set.seed(seed)}

	if (class(dat) == "ExpressionSet"){
	        if (length(grp) == 1) {
			grp <- pData(dat)[, grp]
			dat <- exprs(dat)
		}    
		else{
			dat <- exprs(dat)
		}
	}

	if(is.null(healthy)){
		grp <- make.consecutive.int(grp)
		if (is.null(rownames(dat))) {
			rownames(dat) <- as.character(1:dim(dat)[1])
		}
		if (max(grp) != 1) {
			stop("TSPs can only be calculated for variables with two classes")
		}
		dat <- dat[,c(which(grp==1),which(grp==0))]
		n1 <- length(which(grp==1))
		n2 <- length(which(grp==0))
		grp <- c(rep(1,n1),rep(0,n2))
	}
	else{
		if(!healthy %in% grp){stop("The variable selected for the healthy group does not belong to the dataset")}
		if (length(unique(grp))>2) {
			stop("TSPs can only be calculated for variables with two classes")
		}
		grp[which(grp!=healthy)]<-0
		grp[which(grp==healthy)]<-1
		dat <- dat[,c(which(grp==1),which(grp==0))]
		n1 <- length(which(grp==1))
		n2 <- length(which(grp==0))
		grp <- c(rep(1,n1),rep(0,n2))
	}

	data_boot <-list()
	notin1 <- list()
	notin2<- list()

	for (i in 1:n){	
		random1 <- sample(n1,n1, replace = TRUE)
		random2 <- sample(seq(n1+1,n1+n2),n2, replace = TRUE)
		data_boot[[i]] <- cbind(dat[,random1],dat[,random2])
		notin1[[i]] <- c(seq(1,n1)[!(seq(1,n1) %in% random1)])
		notin2[[i]] <- c(seq(n1+1,n1+n2)[!(seq(n1+1,n1+n2) %in% random2)])
	}

	mean.pred <- function(pred){
		n <- nrow(pred)
		proba <- c()
	
		for(i in 1:n){
			proba[i] <- mean(pred[i,], na.rm=TRUE)
		}
		return(proba)
	}

	prediction <- matrix(nrow=ncol(dat),ncol=n)

	for(m in 1:n){
		cv<- cv(data_boot[[m]], grp, med = med)
		ktsp <- ktspcalc(data_boot[[m]], grp, cv$k, med = med)
		k <- ktsp$k
		pred <- matrix(nrow=ncol(dat), ncol=k)
		for(u in 1:k){
			pred[c(notin1[[m]],notin2[[m]]),u] <- as.numeric(predict(ktsp, dat[,c(notin1[[m]],notin2[[m]])], display=FALSE, select=u))
		}
		prediction[,m] <- mean.pred(pred)
	}

	cutoff <- c(0:10)/10	
	sensitivity <- matrix(nrow=n, ncol=length(cutoff))
	specificity <- matrix(nrow=n, ncol=length(cutoff))

	for(i in 1:n)
		for(cut in cutoff){
			pred <- c()
			pred[prediction[,i]<=cut] <- 0 
			pred[prediction[,i]>cut] <- 1 	
			sensitivity[i,1+10*cut] <- mean(pred[grp==0]==0, na.rm=TRUE)
			specificity[i,1+10*cut] <- mean(pred[grp==1]==1, na.rm=TRUE)
			
	}

	ROC <- list(spec = specificity, sens = sensitivity, n=n)
	class(ROC) <- "ROC"
	return(ROC)
}


