predict.ktsp <- function(object, dat = NULL, select = NULL, display = TRUE,...){
## Compute prediction of the dataset on which the model is based or of a new dataset if one is inserted.

	ktspobj <- object
	grp <- ktspobj$grp
	k <- ktspobj$k
	grplabels <- character(length(grp))
	grplabels[grp==0] <- ktspobj$labels[1]
	grplabels[grp==1] <- ktspobj$labels[2]
	predict <- c()

	if(!is.null(dat) && ktspobj$med == TRUE){
		dat[c(ktspobj$index),] <- dat[c(ktspobj$index),]-ktspobj$med
	}
	
	if(is.vector(dat)){
		dat <- as.matrix(dat,length(dat),1)
	}

	if(!is.null(select) && (select <1 || select > k)){stop("The selected pair of genes is not available.")}

	if(is.null(dat) & !is.null(select)){
		z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
		z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
		table_label <- dimnames(table(z1,grplabels))[[2]]

		p11 <- mean(ktspobj$ktspdat[select,which(grp==0)]<ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
		p21 <- mean(ktspobj$ktspdat[select,which(grp==1)]<ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

		p12 <- mean(ktspobj$ktspdat[select,which(grp==0)]>ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
		p22 <- mean(ktspobj$ktspdat[select,which(grp==1)]>ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

		max <- which.max(c(abs(p11-p21), abs(p12-p22)))

		if(max==1){
			predict[which(z1 == 0)] <- table_label[which.max(table(z1,grplabels)[1,])] 
			predict[which(z1 == 1)] <- table_label[which.max(table(z1,grplabels)[2,])]
		}

		if(max==2){
			predict[which(z2 == 0)] <- table_label[which.max(table(z2,grplabels)[1,])] 
			predict[which(z2 == 1)] <- table_label[which.max(table(z2,grplabels)[2,])]
		}
		return(predict)
	}
	if(is.null(select) && is.null(dat)){

	if(k%%2 == 0){
		k2 <- k-1
		cat("The prediction based on the majority voting procedure cannot be computed for an even value of k. \n")
		cat("The value of k has been reduced to k = ", k2, ". \n")
		k <- k2
	}

	predict <- character(length(grp))
	vote <- numeric(length(grp))
	vote2 <-matrix(nrow=length(grp), ncol=k)
	count <- numeric(length(grp))
	for(i in 1:k){
		z1 <- ktspobj$ktspdat[i,] < ktspobj$ktspdat[(i + k),]
		z2 <- ktspobj$ktspdat[i,] > ktspobj$ktspdat[(i + k),]

		p11 <- mean(ktspobj$ktspdat[i,which(grp==0)]<ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
		p21 <- mean(ktspobj$ktspdat[i,which(grp==1)]<ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

		p12 <- mean(ktspobj$ktspdat[i,which(grp==0)]>ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
		p22 <- mean(ktspobj$ktspdat[i,which(grp==1)]>ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

		max <- which.max(c(abs(p11-p21), abs(p12-p22)))

		if(max==1){
			vote2[which(z1 == 0),i] <- which.max(table(z1,grplabels)[1,]) - 1 
			vote2[which(z1 == 1),i] <- which.max(table(z1,grplabels)[2,]) - 1 
			count[is.na(z1)==FALSE] <- count[is.na(z1)==FALSE]+1
		}

		if(max==2){
			vote2[which(z2 == 0),i] <- which.max(table(z2,grplabels)[1,]) - 1 
			vote2[which(z2 == 1),i] <- which.max(table(z2,grplabels)[2,]) - 1 
			count[is.na(z2)==FALSE] <- count[is.na(z2)==FALSE]+1
		}
	}

	table_label <- labels(table(z1,grplabels))[2]
	vote <- apply(vote2, 1, sum)
	predict[(vote/count < 1/2)] <- table_label$grplabels[1]
	predict[(vote/count > 1/2)] <- table_label$grplabels[2]
	nopred <- which((predict=="")==TRUE)
		if(length(nopred)>0){predict2 <- character(length(nopred))
			if(display==TRUE){
				cat("For the observation(s): ", nopred, " a part of the data is missing \n")
				cat("Their prediction was computed on a subset of the k-TSP\n")
			}
			for(i in 1:length(nopred)){
				vote3 <- numeric(length(nopred))
				vote3 <- vote2[i,is.na(vote2[i,])==FALSE]
				pred <- mean(vote3)
				if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
				if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
				else{
					pred <- mean(vote3[-length(vote3)])
					if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
					if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
				}
			}
			predict[nopred] <- predict2
		}

		return(predict)
	}

	if(!is.null(dat) & !is.null(select)){
		ktspnames <- rownames(ktspobj$ktspdat)[c(select,(select + k))]
		predict <- character(dim(dat)[2])
		vote <- numeric(length(grp))
		vote2 <-matrix(nrow=length(grp), ncol=k)
		count <- numeric(length(grp))
		if(class(dat) == "matrix"){
			if(is.null(rownames(dat))){
				if(display==TRUE){
					cat("No rownames found, using indices \n")}
					ktspnames <- as.numeric(ktspnames)
					if(any(!(ktspnames %in% 1:dim(dat)[1]))){stop("Rownames of new data do not include the TSP names")}
			}
			else{
				if(any(!(ktspnames %in% rownames(dat)))){stop("Rownames of new data do not include the TSP names")}
			}

			z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
			z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
			table_label <- as.vector(dimnames(table(z1,grplabels))[[2]])

			p11 <- mean(ktspobj$ktspdat[select,which(grp==0)]<ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
			p21 <- mean(ktspobj$ktspdat[select,which(grp==1)]<ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

			p12 <- mean(ktspobj$ktspdat[select,which(grp==0)]>ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
			p22 <- mean(ktspobj$ktspdat[select,which(grp==1)]>ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

			max <- which.max(c(abs(p11-p21), abs(p12-p22)))

			if(max==1){
				w <- dat[ktspnames[1],] < dat[ktspnames[2],]
				predict[which(w == 0)] <- table_label[which.max(table(z1,grplabels)[1,])]
				predict[which(w == 1)] <- table_label[which.max(table(z1,grplabels)[2,])]
			}

			if(max==2){
				w <- dat[ktspnames[1],] > dat[ktspnames[2],]
				predict[which(w == 0)] <- table_label[which.max(table(z2,grplabels)[1,])] 
				predict[which(w == 1)] <- table_label[which.max(table(z2,grplabels)[2,])]
			}
			nopred <- which((predict=="")==TRUE)

			if(length(nopred)>0){
				if(display==TRUE){
					cat("For the observation(s): ", nopred, " a part of the data is missing \n")
					cat("Their prediction for the selected pair is not possible\n")}
				}
				return(predict)
			}
			if(class(dat) == "ExpressionSet"){
				if(is.null(featureNames(dat))){
					if(display==TRUE){
						cat("No featureNames info found, using indices \n")
					}
					ktspnames <- as.numeric(ktspnames)
					if(any(!(ktspnames %in% 1:dim(exprs(dat))[1]))){stop("Rownames of new data do not include the TSP names")}
					dat <- exprs(dat)
				}
			else{
				if(any(!(ktspnames %in% featureNames(dat)))){stop("Rownames of new data do not include the TSP names")}
					genenames <- featureNames(dat)
					dat <- exprs(dat)
					rownames(dat) <- genenames
			}
			z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
			z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
			table_label <- dimnames(table(z1,grplabels))[[2]]

			z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
			z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
			table_label <- as.vector(dimnames(table(z1,grplabels))[[2]])

			p11 <- mean(ktspobj$ktspdat[select,which(grp==0)]<ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
			p21 <- mean(ktspobj$ktspdat[select,which(grp==1)]<ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

			p12 <- mean(ktspobj$ktspdat[select,which(grp==0)]>ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
			p22 <- mean(ktspobj$ktspdat[select,which(grp==1)]>ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

			max <- which.max(c(abs(p11-p21), abs(p12-p22)))

			if(max==1){
				w <- dat[ktspnames[1],] < dat[ktspnames[2],]
				predict[which(w == 0)] <- table_label[which.max(table(z1,grplabels)[1,])]
				predict[which(w == 1)] <- table_label[which.max(table(z1,grplabels)[2,])]
			}

			if(max==2){
				w <- dat[ktspnames[1],] > dat[ktspnames[2],]
				predict[which(w == 0)] <- table_label[which.max(table(z2,grplabels)[1,])] 
				predict[which(w == 1)] <- table_label[which.max(table(z2,grplabels)[2,])]
			}
			nopred <- which((predict=="")==TRUE)
			if(length(nopred)>0){
				if(display==TRUE){
					cat("For the observation(s): ", nopred, " a part of the data is missing \n")
					cat("Their prediction for the selected pair is not possible\n")
				}
			}
			return(predict)
		}
	}

	else{

		if(k%%2 == 0){
			k2 <- k-1
			if(display==TRUE){
				cat("The prediction based on the k-TSP cannot be computed for an even value of k. \n")
				cat("The value of k has been reduced to k = ", k2, ". \n")
			}
			k <- k2
		}
		ktspnames <- rownames(ktspobj$ktspdat)
		if(class(dat) == "matrix"){
			if(is.null(rownames(dat))){
				if(display==TRUE){
					cat("No rownames found, using indices \n")
				}
				ktspnames <- as.numeric(ktspnames)
				if(any(!(ktspnames %in% 1:dim(dat)[1]))){stop("Rownames of new data do not include the TSP names")}
			}

			else{
				if(any(!(ktspnames %in% rownames(dat)))){stop("Rownames of new data do not include the TSP names")}
			}
			predict <- character(dim(dat)[2])
			vote <- numeric(dim(dat)[2])
			vote2 <-matrix(nrow=dim(dat)[2], ncol=k)
			count <- numeric(dim(dat)[2])
			for(i in 1:k){
				z1 <- ktspobj$ktspdat[i,] < ktspobj$ktspdat[(i + k),]
				z2 <- ktspobj$ktspdat[i,] > ktspobj$ktspdat[(i + k),]

				p11 <- mean(ktspobj$ktspdat[i,which(grp==0)]<ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
				p21 <- mean(ktspobj$ktspdat[i,which(grp==1)]<ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

				p12 <- mean(ktspobj$ktspdat[i,which(grp==0)]>ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
				p22 <- mean(ktspobj$ktspdat[i,which(grp==1)]>ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

				max <- which.max(c(abs(p11-p21), abs(p12-p22)))

				if(max==1){
					w <- dat[ktspnames[i],] < dat[ktspnames[i+k],]
					vote2[which(w == 0),i] <- which.max(table(z1,grplabels)[1,]) - 1 
					vote2[which(w == 1),i] <- which.max(table(z1,grplabels)[2,]) - 1 
					count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
				}

				if(max==2){
					w <- dat[ktspnames[i],] > dat[ktspnames[i+k],]
					vote2[which(w == 0),i] <- which.max(table(z2,grplabels)[1,]) - 1 
					vote2[which(w == 1),i] <- which.max(table(z2,grplabels)[2,]) - 1 
					count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
				}
			}
			table_label <- labels(table(z1,grplabels))[2]
			vote <- apply(vote2, 1, sum)
			predict[(vote/count < 1/2)] <- table_label$grplabels[1]
			predict[(vote/count > 1/2)] <- table_label$grplabels[2]
			nopred <- which((predict=="")==TRUE)

#If predictions were not possible for some observations, the prediction will be based on a restricted list of the pairs in the k-TSP.

			if(length(nopred)>0){
				predict2 <- character(length(nopred))
				if(display==TRUE){
					cat("For the observation(s): ", nopred, " a part of the data is missing \n")
					cat("Their prediction was computed on a subset of the k-TSP\n")
				}
			for(i in 1:length(nopred)){
				vote3 <- numeric(length(nopred))
				vote3 <- vote2[nopred[i],is.na(vote2[nopred[i],])==FALSE]
				if(length(vote3)<1){if(display==TRUE){cat("Prediction not possible for observation", nopred[i], " \n")}}
				else{
					pred <- mean(vote3)
					if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
					if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
					else{
						if(length(vote3)<2){if(display==TRUE){cat("Prediction not possible for observation", nopred[i], "\n")}}
						else{
							pred <- mean(vote3[-length(vote3)])
							if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
							if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
						}
					}
				}
			}
			predict[nopred] <- predict2
			}
			return(predict)
		}
	
		if(class(dat) == "ExpressionSet"){
			if(is.null(featureNames(dat))){
				if(display==TRUE){
					cat("No featureNames info found, using indices \n")}
				ktspnames <- as.numeric(ktspnames)
				if(any(!(ktspnames %in% 1:dim(exprs(dat))[1]))){stop("Rownames of new data do not include the TSP names")}
				dat <- exprs(dat)
			}

			else{
			if(any(!(ktspnames %in% featureNames(dat)))){stop("Rownames of new data do not include the TSP names")}
				genenames <- featureNames(dat)
				dat <- exprs(dat)
				rownames(dat) <- genenames
			}
			predict <- character(dim(dat)[2])
			vote <- numeric(dim(dat)[2])
			vote2 <-matrix(nrow=dim(dat)[2], ncol=k)
			count <- numeric(dim(dat)[2])
			for(i in 1:k){
				z1 <- ktspobj$ktspdat[i,] < ktspobj$ktspdat[(i + k),]
				z2 <- ktspobj$ktspdat[i,] > ktspobj$ktspdat[(i + k),]

				p11 <- mean(ktspobj$ktspdat[i,which(grp==0)]<ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
				p21 <- mean(ktspobj$ktspdat[i,which(grp==1)]<ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

				p12 <- mean(ktspobj$ktspdat[i,which(grp==0)]>ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
				p22 <- mean(ktspobj$ktspdat[i,which(grp==1)]>ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

				max <- which.max(c(abs(p11-p21), abs(p12-p22)))

				if(max==1){
					w <- dat[ktspnames[i],] < dat[ktspnames[i+k],]
					vote2[which(w == 0),i] <- which.max(table(z1,grplabels)[1,]) - 1 
					vote2[which(w == 1),i] <- which.max(table(z1,grplabels)[2,]) - 1 
					count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
				}

				if(max==2){
					w <- dat[ktspnames[i],] > dat[ktspnames[i+k],]
					vote2[which(w == 0),i] <- which.max(table(z2,grplabels)[1,]) - 1 
					vote2[which(w == 1),i] <- which.max(table(z2,grplabels)[2,]) - 1 
					count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
				}
			}
			table_label <- labels(table(z1,grplabels))[2]
			vote <- apply(vote2, 1, sum)
			predict[(vote/count < 1/2)] <- table_label$grplabels[1]
			predict[(vote/count > 1/2)] <- table_label$grplabels[2]
			nopred <- which((predict=="")==TRUE)
			if(length(nopred)>0){
				predict2 <- character(length(nopred))
				if(display==TRUE){
					cat("For the observation(s): ", nopred, " a part of the data is missing \n")
					cat("Their prediction was computed on a subset of the k-TSP\n")
				}
				for(i in 1:length(nopred)){
					vote3 <- numeric(length(nopred))
					vote3 <- vote2[nopred[i],is.na(vote2[nopred[i],])==FALSE]
					if((length(vote3)<1) & (display == TRUE)){cat("Prediction not possible for observation", nopred[i], " \n")}
					else{
						pred <- mean(vote3)
						if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
						if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
						else{
							if((length(vote3)<2) & (display == TRUE)){cat("Prediction not possible for observation", nopred[i], "\n")}
							else{
								pred <- mean(vote3[-length(vote3)])
								if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
								if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
							}
						}
					}
				}
				predict[nopred] <- predict2
			}
			return(predict)
		}	 
	}
}

