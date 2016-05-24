seqplot.rf <- function(seqdata, k=floor(nrow(seqdata)/10), diss, sortv=NULL,
						ylab=NA, yaxis=FALSE, title=NULL, ...){
	
	return(seqplot.rf_internal(seqdata, k=k, diss=diss, sortv=sortv,
						ylab=ylab, yaxis=yaxis, main=title, ...))
}
seqplot.rf_internal <- function(seqdata, k=floor(nrow(seqdata)/10), diss, sortv=NULL,
						use.hclust=FALSE, hclust_method="ward.D", use.quantile=FALSE,
						yaxis=FALSE, main=NULL, ...){
	
	message(" [>] Using k=", k, " frequency groups")
	
	
	#Extract medoid, possibly weighted
	gmedoid.index <- disscenter(diss, medoids.index="first")
	
	gmedoid.dist <-diss[, gmedoid.index] #Extract distance to general medoid

	##Vector where distance to k medoid will be stored
	kmedoid.dist <- rep(0, nrow(seqdata))
	#index of the k-medoid for each sequence
	kmedoid.index <- rep(0, nrow(seqdata))
	#calculate qij - distance to frequency group specific medoid within frequency group
	if(is.null(sortv) && !use.hclust){
		sortv <- cmdscale(diss, k = 1)
	
	}
	if(!is.null(sortv)){
		ng <- nrow(seqdata) %/% k
		r <- nrow(seqdata) %% k
		n.per.group <- rep(ng, k)
		if(r>0){
			n.per.group[order(runif(r))] <- ng+1
		}
		mdsk <- rep(1:k, n.per.group)
		mdsk <- mdsk[rank(sortv, ties.method = "random")]
	}else{
		hh <- hclust(as.dist(diss), method=hclust_method)
		mdsk <- factor(cutree(hh, k))
		medoids <- disscenter(diss, group=mdsk, medoids.index="first")
		medoids <- medoids[levels(mdsk)]
		#ww <- xtabs(~mdsk)
		mds <- cmdscale(diss[medoids, medoids], k=1)
		mdsk <- as.integer(factor(mdsk, levels=levels(mdsk)[order(mds)]))
	}
	kun <- length(unique(mdsk))
	if(kun!=k){
		warning(" [>] k value was adjusted to ", kun)
		k <- kun
		mdsk <- as.integer(factor(mdsk, levels=sort(unique(mdsk))))
	}
	#sortmds.seqdata$mdsk<-c(rep(1:m, each=r+1),rep({m+1}:k, each=r))
	##pmdse <- 1:k
	#pmdse20<-1:20
	
	##for each k
	for(i in 1:k){
		##Which individuals are in the k group
		ind <- which(mdsk==i)
		if(length(ind)==1){
			kmedoid.dist[ind] <- 0
			##Index of the medoid sequence for each seq
			kmedoid.index[ind] <- ind
		}else{
			dd <- diss[ind, ind]
			##Indentify medoid
			kmed <- disscenter(dd, medoids.index="first")
			##Distance to medoid for each seq
			kmedoid.dist[ind] <- dd[, kmed]
			##Index of the medoid sequence for each seq
			kmedoid.index[ind] <- ind[kmed]
		}
		##Distance matrix for this group
		
	}

	##Attribute to each sequences the medoid sequences
	seqtoplot <- seqdata[kmedoid.index, ]
	
	##Correct weights to their original weights (otherwise we use the medoid weights)
	attr(seqtoplot, "weights") <- NULL
	opar <- par(mfrow=c(1,2), oma=c(3,0,(!is.null(main))*3,0), mar=c(1, 1, 2, 0))
	on.exit(par(opar))
	seqIplot(seqtoplot, withlegend=FALSE, sortv=mdsk, yaxis=yaxis, title="Sequences medoids", ...)
	##seqIplot(seqtoplot, withlegend=FALSE, sortv=mdsk)
	heights <- xtabs(~mdsk)/nrow(seqdata)
	at <- (cumsum(heights)-heights/2)/sum(heights)*length(heights)
	if(!yaxis){
		par(yaxt="n")
	}
	
	boxplot(kmedoid.dist~mdsk, horizontal=TRUE, width=heights, frame=FALSE,  main="Dissimilarities to medoid", ylim=range(as.vector(diss)), at=at)
	
	#calculate R2
	R2 <-1-sum(kmedoid.dist^2)/sum(gmedoid.dist^2)
	#om K=66 0.5823693
	
	
	#calculate F
	ESD <-R2/(k-1) # averaged explained variance
	USD <-(1-R2)/(nrow(seqdata)-k) # averaged explained variance
	Fstat <- ESD/USD
	
	message(" [>] Pseudo/median-based-R2: ", format(R2))
	message(" [>] Pseudo/median-based-F statistic: ", format(Fstat))
	##cat(sprintf("Representation quality: R2=%0.2f F=%0.2f", R2, Fstat))
	title(main=main, outer=TRUE)
	title(sub=sprintf("Representation quality: R2=%0.2f and F=%0.2f", R2, Fstat), outer=TRUE, line=2)
}
