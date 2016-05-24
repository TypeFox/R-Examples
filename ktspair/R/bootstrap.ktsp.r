

bootstrap.ktsp <- function(dat, grp, k = NULL, seed = NULL, n = 500, length = 40, display = FALSE, med = FALSE){

	if(!is.null(seed)){set.seed(seed)}

	if (class(dat) == "ExpressionSet"){
        	genenames <- as.character(1:dim(exprs(dat))[1])
     		if (!is.null(featureNames(dat))) {
				genenames <- featureNames(dat)
        	}
        	if (length(grp) == 1) {
			grp <- make.consecutive.int(pData(dat)[, grp])
			dat <- exprs(dat)
			rownames(dat) <- genenames
				if (max(grp) != 1) {
					stop("TSPs can only be calculated for variables with two classes")
				}
		}
		else{
			grp <- make.consecutive.int(grp)
			dat <- exprs(dat)
			rownames(dat) <- genenames
			if (max(grp) != 1) {
				stop("TSPs can only be calculated for variables with two classes")
			}
		}
	}
	else {
        	grp <- make.consecutive.int(grp)
		genenames <- as.character(1:dim(dat)[1])
        	if (!is.null(rownames(dat))) {
			genenames <- rownames(dat)
        	}
        	if (max(grp) != 1) {
			stop("TSPs can only be calculated for variables with two classes")
        	}
	}
			

	dat <- dat[,c(which(grp==0),which(grp==1))]
	n1 <- length(which(grp==0))
	n2 <- length(which(grp==1))

	data_boot1 <-list()

	for (i in 1:n){	
		sample <-cbind(dat[,sample(n1,n1, replace = TRUE)],dat[,sample(seq(n1+1,n1+n2),n2, replace = TRUE)])
		data_boot1[[i]] <- sample
	}

	cv <- cv(dat, grp, length = length, display = display, med = med)
	ktsp <- ktspcalc(dat, grp, cv$k, display = display, med = med)

	score <- matrix(nrow=n,ncol=9)
	index <- matrix(nrow=n,ncol=18)
	k_value <- c()

	if(is.null(k)){
		for(i in 1:n){
			ktsp<-ktspcalc(data_boot1[[i]], grp, length= length, display = display, med = med)
			score[i,1:length(ktsp$ktspscore)] <- ktsp$ktspscore
			index[i,1:(2*ktsp$k)] <- as.vector(t(ktsp$index)) 
			k_value[i] <- ktsp$k
		}
	}
	else {
		for(i in 1:n){
			ktsp<-ktspcalc(data_boot1[[i]], grp, k, length = length, display = display, med = med)
			score[i,1:length(ktsp$ktspscore)] <- ktsp$ktspscore
			index[i,1:(2*ktsp$k)] <- as.vector(t(ktsp$index)) 
			k_value[i] <- ktsp$k
		}
	}
	bootstrap <- list(score = score, index = index, k_value = k_value, k = k, n = n, genenames = genenames, ktsp = ktsp)
	class(bootstrap) <- "bootstrap"
	return(bootstrap)
}


