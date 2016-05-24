

ktspcalc2 <- function(indice, dat, grp, length = 40, cross =5, display = display, med = FALSE, performance = FALSE, seed = NULL, healthy = NULL){

	accuracy <- NULL
	accuracy_k <- NULL
	sensitivity <- NULL
	specificity <- NULL

	if(!is.null(healthy)){
		if(!(healthy %in% grp)){stop("The group called ",healthy," is not present in the variable grp \n")}
		else performance <- TRUE
	}


	if (class(dat) == "ExpressionSet"){
        	genenames <- as.character(1:dim(exprs(dat))[1])
     		if (!is.null(featureNames(dat))) {
				genenames <- featureNames(dat)
        	}
        	if (length(grp) == 1) {
			labels <- as.character(unique(pData(dat)[, grp])[order(unique(grp))])
			if(!is.null(healthy)){healthy.position <- which((pData(dat)[, grp])==healthy)[1]}
			grp <- make.consecutive.int(pData(dat)[, grp])
			if(!is.null(healthy)){healthy <- grp[healthy.position]}
			dat <- exprs(dat)
			rownames(dat) <- genenames
				if (max(grp) != 1) {
					stop("TSPs can only be calculated for variables with two classes")
				}
		}
		else{
			labels <- as.character(unique(grp)[order(unique(grp))])
			if(!is.null(healthy)){healthy.position <- which(grp==healthy)[1]}
			grp <- make.consecutive.int(grp)
			if(!is.null(healthy)){healthy <- grp[healthy.position]}
			dat <- exprs(dat)
			rownames(dat) <- genenames
			if (max(grp) != 1) {
				stop("TSPs can only be calculated for variables with two classes")
			}
		}
	}
	else {
	        labels <- as.character(unique(grp)[order(unique(grp))])
		if(!is.null(healthy)){healthy.position <- which(grp==healthy)[1]}
		grp <- make.consecutive.int(grp)
		if(!is.null(healthy)){healthy <- grp[healthy.position]}
		genenames <- as.character(1:dim(dat)[1])
		rownames(dat) <- genenames
        	if (!is.null(rownames(dat))) {
			genenames <- rownames(dat)
        	}
        	if (max(grp) != 1) {
			stop("TSPs can only be calculated for variables with two classes")
        	}
	}
	
	if(med){
		median_0 <- c(apply(dat[,which(grp==1)], 1, function(x) median(x,na.rm=TRUE)))
		median_1 <- c(apply(dat[,which(grp==0)], 1, function(x) median(x,na.rm=TRUE)))
		mean_med <- (median_0+median_1)/2
		dat <- dat-mean_med
	}
			
	if(performance){
		cv <- cv2(indice=indice,dat, grp, cross = cross, display = display, length = length, healthy = healthy, seed = seed)
		accuracy <- cv$accuracy
		if(!is.null(healthy)){
			sensitivity <- cv$sensitivity
			specificity <- cv$specificity
		}
	}

	k <- length(indice)/2
	index <- matrix(indice,ncol=2, nrow=k, byrow=TRUE) 

	p11 <- apply(matrix(c(dat[c(index[,1]),which(grp==0)]<dat[c(index[,2]),which(grp==0)]),k,length(which(grp==0)), byrow=F), 1, function(x) mean(na.omit(x)))
	p21 <- apply(matrix(c(dat[c(index[,1]),which(grp==1)]<dat[c(index[,2]),which(grp==1)]),k,length(which(grp==1)), byrow=F), 1, function(x) mean(na.omit(x)))

	p12 <- apply(matrix(c(dat[c(index[,1]),which(grp==0)]>dat[c(index[,2]),which(grp==0)]),k,length(which(grp==0)), byrow=F), 1, function(x) mean(na.omit(x)))
	p22 <- apply(matrix(c(dat[c(index[,1]),which(grp==1)]>dat[c(index[,2]),which(grp==1)]),k,length(which(grp==1)), byrow=F), 1, function(x) mean(na.omit(x)))

	score1 <- c(abs(p11-p21))
	score2 <- c(abs(p12-p22))

	score <- rbind(score1, score2)
	ktspscore <- apply(score, 2, max)

	rank <- rank_na(dat, NA)
	
	r1 <- apply(matrix(c(rank[c(index[,1]),which(grp==0)]-rank[c(index[,2]),which(grp==0)]),k,length(which(grp==0)), byrow=T), 1, function(x) mean(na.omit(x)))
	r2 <- apply(matrix(c(rank[c(index[,1]),which(grp==1)]-rank[c(index[,2]),which(grp==1)]),k,length(which(grp==1)), byrow=T), 1, function(x) mean(na.omit(x)))
	
	rankscore <- c(abs(r1-r2))

	order_tsp <- ordertsp(ktspscore, rankscore)

	ktspdat <- dat[c(index[,1][order_tsp],index[,2][order_tsp]),]


	index <- matrix(c(index[order_tsp,]),ncol=2, nrow=k, byrow=FALSE)

	if(med) med <- mean_med[c(index)]

	ktsp <- list(index = index, ktspscore = ktspscore[order_tsp], grp = grp, ktspdat = ktspdat, k = k, labels =labels, rankscore = rankscore[order_tsp],  accuracy =accuracy, sensitivity = sensitivity, specificity = specificity, med = med)
	class(ktsp) <- "ktsp"
	return(ktsp)

}


