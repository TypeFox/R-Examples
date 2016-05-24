makeDyadic <-
function(x,directed=FALSE, show.progress=5){
	
	# we take out all years with only one country since there is no dyad there
	take.out <- which(table(dataOrig[,2])==1)
	take.out.n <- as.numeric(names(take.out))
	rows.t.o <- which(dataOrig[,2] %in% take.out.n)
	ifelse(length(rows.t.o)==0,dataOrig <- dataOrig,dataOrig <- dataOrig[-rows.t.o,])


	# First step is to slice up the for each year
	# (second variable has to be year)	
	n.time <- length(table(dataOrig[,2]))
	Years <- c(table(dataOrig[,2]))
	years <- as.numeric(names(Years))
	dataList <- list()
	for (j in years){
		aa <- paste("data",j,sep="")
		dataList[[aa]] <- dataOrig[dataOrig[,2]==j,]
	}


	# After having created a list in which each element has all the data from a given year, we
	# go and take each of these elements and match every year with every other one (note: we
	# have directed dyads)



	dataListdyad <- list()

	if (directed==1)	for (w in 1:n.time){
			name.dyadYEAR <- paste("data.",w,".year",sep="")
			N1 <- dim(dataList[[w]])[1]
			N2 <- dim(dataList[[w]])[2]
			dataListdyad[[name.dyadYEAR]] <- data.frame(matrix(NA,N1*(N1-1),2*N2))
			Veci <- c(1:N1)
			counter <- 0
			for (k in 1:N1){
				veci <- Veci[Veci!=k]
				for (i in veci){
					counter <- counter+1
					dataListdyad[[w]][counter,] <- c(dataList[[w]][k,],dataList[[w]][i,])	
				}
			}
		if (w %% show.progress == 0) writeLines(paste(w,"out of",n.time, "years have been transformed - stay calm"))
		}
	
	if (directed==0)	for (w in 1:n.time){
			name.dyadYEAR <- paste("data.",w,".year",sep="")
			N1 <- dim(dataList[[w]])[1]
			N2 <- dim(dataList[[w]])[2]
			dataListdyad[[name.dyadYEAR]] <- data.frame(matrix(NA,N1*(N1-1)/2,2*N2))
			Veci <- c(1:N1)
			counter <- 0
			for (k in 1:N1){
				veci <- Veci[Veci!=k&Veci>k]
				for (i in veci){
					counter <- counter+1
					dataListdyad[[w]][counter,] <- c(dataList[[w]][k,],dataList[[w]][i,])	
				}
			}
		if (w %% show.progress == 0) writeLines(paste(w,"out of",n.time, "years have been transformed - stay calm"))
		}

	# Now, you find your dyadic data set as a list where each element of the list has all the 
	# dyads from a specific year. In a final step we bind all these elements together to have 
	# one big matrix and add correct variable names

	BigData <- dataListdyad[[1]]
		for (k in 2:n.time){
			BigData <- rbind(BigData,dataListdyad[[k]])
		if (k %% show.progress == 0) writeLines(paste(k,"years already combined"))
		}

	name1 <- colnames(dataOrig)
	name.11 <- rep(NA,length(name1))
	for(i in 1:length(name1)){
		name.11[i] <-  paste("sen",name1[i],sep="_")
	}

	name1 <- colnames(dataOrig)
	name.12 <- rep(NA,length(name1))
	for(i in 1:length(name1)){
		name.12[i] <-  paste("rec",name1[i],sep="_")
	}

	colnames(BigData) <- c(name.11, name.12)
	
	return(BigData)
}
