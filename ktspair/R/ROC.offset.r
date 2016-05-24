

ROC.offset <- function(dat, grp, n = 200, healthy = NULL, seed = NULL, para1 = 200, para2 = 1, mult.cutoff = FALSE, length = 40, display = FALSE, med = FALSE){

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

	mean.pred <- function(single_pred){
		n <- nrow(single_pred[[1]])
		m <- ncol(single_pred[[1]])
		resp <- matrix(ncol=m, nrow=n)
	
		for(i in 1:n){
			for(j in 1:m){
				vect_mean <- c()
				for(l in 1:length(single_pred)){
					vect_mean[l] <- single_pred[[l]][i,j]
				}
				resp[i,j] <- mean(vect_mean, na.rm=TRUE)
			}
		}
		return(resp)
	}

	if(mult.cutoff==FALSE){	

		spec_list <- matrix(nrow=n,ncol=2*(para1+para2)+1)
		sens_list <- matrix(nrow=n,ncol=2*(para1+para2)+1)

		for(m in 1:n){
			cv<- cv(data_boot[[m]], grp, med = med, length = length, display = display)

			ktsp <- ktspcalc(data_boot[[m]], grp, cv$k, med = med, length = length, display = display)
			k <- ktsp$k

			dat <- data_boot[[1]]
			data <- dat
			i <- ktsp$index[1,1]
			j <- ktsp$index[1,2]

			maxk <- c()
			for(i in 1:k){
				maxk[i] <- max(abs(dat[ktsp$index[i,1],]-dat[ktsp$index[i,2],]), na.rm=TRUE)
			}

			A <- max(maxk)
			D <- A/para1
			seq <- seq(-A-para2*D,A+para2*D,D)
	
			prediction <- c(length(grp))
			spec <-c()
			sens <- c()
			sign <- c()

			for(i in 1:k){
				p1 <- mean(data_boot[[m]][ktsp$index[i,1],which(grp==0)]<data_boot[[m]][ktsp$index[i,2],which(grp==0)], na.rm=TRUE)
				p2 <- mean(data_boot[[m]][ktsp$index[i,1],which(grp==1)]<data_boot[[m]][ktsp$index[i,2],which(grp==1)], na.rm=TRUE)
				if(p1>p2){sign[i]<-1}
				else{sign[i]<--1}
			}

			for(i in 1:length(seq)){
				data[c(ktsp$index[,1]),] <- dat[c(ktsp$index[,1]),] + sign * seq[i]
				prediction[c(notin1[[m]],notin2[[m]])] <- as.numeric(predict(ktsp, data[,c(notin1[[m]],notin2[[m]])], display=FALSE))
				spec[i] <- mean((prediction[1:n1]==grp[1:n1]), na.rm=TRUE)
				sens[i] <- mean((prediction[(n1+1):(n1+n2)]==grp[(n1+1):(n1+n2)]), na.rm=TRUE)
			}

			spec_list[m,] <- spec
			sens_list[m,] <- sens
		}

	}

	else{
		cutoff <- c(0.25, 0.5, 0.75)
		sens_list <- list(length=3)
		spec_list <- list(length=3)

		for(i in 1:3){
			spec_list[[i]] <- matrix(nrow=n,ncol=2*(para1+para2)+1)
			sens_list[[i]] <- matrix(nrow=n,ncol=2*(para1+para2)+1)
		}


		for(m in 1:n){
			cv<- cv(data_boot[[m]], grp, length = length, display = display, med = med)
			ktsp <- ktspcalc(data_boot[[m]], grp, cv$k, length = length, display = display, med = med)
			k <- ktsp$k

			data <- dat
			i <- ktsp$index[1,1]
			j <- ktsp$index[1,2]
			maxk <- c()
			for(i in 1:k){
				maxk[i] <- max(abs(dat[ktsp$index[i,1],]-dat[ktsp$index[i,2],]), na.rm=TRUE)
			}
			A <- max(maxk)
			D <- A/para1
			seq <- seq(-A-para2*D,A+para2*D,D)
	
			single_pred <- list(length=ktsp$k)

			for(i in 1:k){
				single_pred[[i]] <- matrix(nrow=length(seq), ncol=length(grp))
			}
	
			for(l in 1:length(seq)){
				data[c(ktsp$index[,1]),] <- dat[c(ktsp$index[,1]),] + seq[l]		
				for(u in 1:k){
					single_pred[[u]][l,c(notin1[[m]],notin2[[m]])] <- as.numeric(predict(ktsp, data[,c(notin1[[m]],notin2[[m]])],select = u, display = display))
				}
			}	

			predict <- matrix(nrow=length(seq), ncol=length(grp))
			predict <- mean.pred(single_pred)


			for(M in cutoff){
				sens <- c(length=length(seq))
				spec <- c(length=length(seq))
				
				for(i in 1:length(seq)){
					prediction <- c(length=length(grp))
					prediction[which(predict[i,]<M)] <- rep(0,length(which(predict[i,]<M)))
					prediction[which(predict[i,]>=M)] <- rep(1,length(which(predict[i,]>=M)))
		
					spec[i] <- mean((prediction[1:n1]==grp[1:n1]), na.rm=TRUE)
					sens[i] <- mean((prediction[(n1+1):(n1+n2)]==grp[(n1+1):(n1+n2)]), na.rm=TRUE)
				}
				spec_list[[4*M]][m,] <- spec
				sens_list[[4*M]][m,] <- sens
			}
		}
	}
	ROC <- list(spec = spec_list, sens = sens_list, mult.cutoff = mult.cutoff, n=n)
	class(ROC) <- "ROC"
	return(ROC)
}

