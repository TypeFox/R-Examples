
bootstrap.graphic.ktsp <- function(bootstrap, para1 = 0, para2 = 0, title = NULL, mtext = NULL){

	index <- bootstrap$index
	k_value <- bootstrap$k_value
	k <- bootstrap$k
	n <- bootstrap$n
	genenames <- bootstrap$genenames
	ktsp <- bootstrap$ktsp

	if(is.null(k)){
		layout(rbind(c(1,1,2,2),c(4,3,3,4)))
		par(oma=c(0,0,3,0))
		par(mar=c(4,4,7,2))
	}

	else {
		layout(c(1,1,2,2))
		par(oma=c(0,0,3,0))
		par(mar=c(4,4,7,2))
	}

	unique <- unique(as.vector(index))
	count <- c(length=length(unique))

	for(i in 1:length(unique)){
		count[i] <- length(which(index == unique[i]))
	}

	count1 <- count/n

	plot(unique, count1 , xlab="Genes", ylab="Frequency", type = "h", ylim=c(0,max(count1)*1.25))

	title("Frequency of single genes in the k-tsp", line=1)
	text(unique[count1>para1], count1[count1>para1], labels = genenames[unique[count1>para1]], pos=3)


## Graph of the pairs that appear the most often in the k-tsp

	index2 <- matrix(nrow=n, ncol=9)
	names2 <- matrix(nrow=n, ncol=9)

	for(i in 1:dim(index)[1]){
		for(j in 1:9){
			index2[i,j] <- paste (index[i,2*j-1],index[i,2*j], sep="/")
			names2[i,j] <- paste (genenames[index[i,2*j-1]],genenames[index[i,2*j]], sep="/")
			if(index2[i,j]=="NA/NA"){index2[i,j] <-NA}
			if(names2[i,j]=="NA/NA"){names2[i,j] <-NA}
		}
	}


	index2 <- index2[!is.na(index2)]
	names2 <- names2[!is.na(names2)]
	unique <- unique(as.vector(index2))
	names <- unique(as.vector(names2))
	count <- c(length=length(unique))

	for(i in 1:length(unique)){
		count[i] <- length(which(index2 == unique[i]))/n
	}


	count1 <- count[count>para2]
	unique2 <- unique[count>para2]
	names <- names[count>para2]

	indice <- matrix(nrow=dim(ktsp$index)[1],ncol=dim(ktsp$index)[2])
	for(i in 1:dim(ktsp$index)[1]){
		if(paste(ktsp$index[i,1], ktsp$index[i,2], sep="/") %in% unique2){
			indice[i,] <- which(unique2 == paste(ktsp$index[i,1], ktsp$index[i,2], sep="/"))
		}
	}


	pos <- rep(c(1,3),length(unique2)/2+1)

	plot(1:length(unique2), count1 , xlab="Pair of genes",  ylab="Frequency", type = "p", ylim=c(0,max(count1)*1.25), xlim=c(0,length(unique2)*1.1), pch=20)

	title("Frequency of gene pairs in the k-tsp", line=1)
	text(1:length(unique2), count1, labels = names, pos=pos)

	for(j in 1:length(indice)){
		if(is.na(indice[j])==FALSE){
			points(indice[j], count1[indice[j]],pch=20, col="red")
		}
	}

	if(is.null(k)){## Histogram for the values of k
		hist(k_value, breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5), freq=FALSE, xlim=c(0,9), ylim=c(0,1),
		xlab="k", axes=FALSE, main=NA, col = "blue", ylab="Frequency")
		axis(side = 1, at = c(1,3,5,7,9), labels = c(1,3,5,7,9))
		axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
		title("Histogram of the values of k in the k-TSP", line=1)
	}

	if(is.null(title)){
	title("Summary of the results for the k-TSP", line=-1, outer=TRUE, cex.main=2)
	}

	else{
	title(title, line=-1.5, outer=TRUE, cex.main=2)
	}
	if(!is.null(mtext)){
		mtext(mtext, side=3, line=-4, outer=TRUE)
	}
}


