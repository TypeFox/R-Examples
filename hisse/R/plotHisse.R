
######################################################################################################################################
######################################################################################################################################
### Plotting function for HiSSE
######################################################################################################################################
######################################################################################################################################

plot.hisse.states <- function(x, rate.param, do.observed.only=TRUE, rate.colors=NULL, state.colors=NULL, edge.width.rate=5, edge.width.state=2, type="fan", rate.range=NULL, show.tip.label=TRUE, fsize=1.0, lims.percentage.correction=0.001, legend="tips", legend.position=c(0, 0.2, 0, 0.2), legend.cex=0.4, legend.kernel.rates="auto", legend.kernel.states="auto", legend.bg="cornsilk3", ...) {
	hisse.results <- x
	if(class(hisse.results)=="hisse.states") { #we have to make a list so we can run this generally
		if(is.null(hisse.results$aic)){
			#If a user forgot to include the aic, then we add a random value in for them
			hisse.results$aic = 42
		}
		tmp.list <- list()
		tmp.list[[1]] <- hisse.results
		hisse.results <- tmp.list	
	}
	par(fig=c(0,1, 0, 1), new=FALSE)
	if(is.null(rate.colors)) {
		rate.colors <- c("blue", "red")	
	}
	if(is.null(state.colors)) {
		state.colors <- c("white", "black")
	}
	rates.tips <- ConvertManyToRate(hisse.results, rate.param, "tip.mat")
	rates.internal <- ConvertManyToRate(hisse.results, rate.param, "node.mat")
	states.tips <- NA
	states.internal <- NA
	if (do.observed.only) {
		states.tips <- ConvertManyToBinaryState(hisse.results, "tip.mat")
		states.internal <- ConvertManyToBinaryState(hisse.results, "node.mat")
	} else {
		stop("So far we can easily plot just the binary observed state; if you want to plot the hidden states, use a different function")
	}
	tree.to.plot <- hisse.results[[1]]$phy
	if(!show.tip.label) { #this is b/c you cannot suppress plotting tip labels in phytools plotSimmap
		rep.rev <- function(x, y) {
			result<-paste(rep(y,x), collapse="", sep="")
			return(result)	
		}
		tree.to.plot$tip.label <- sapply(sequence(length(tree.to.plot$tip.label)), rep.rev, " ")
		fsize=0
	}
#	rate.tree <- contMapGivenAnc(tree=hisse.object$phy, x=ConvertToRate(hisse.object$tip.mat, rate.vector= rate.vector), plot=FALSE, anc.states=ConvertToRate(hisse.object$node.mat, rate.vector= rate.vector), ...)
	rate.lims <- range(c(rates.tips, rates.internal))
	if(!is.null(rate.range)) {
		if(min(rate.range) > min(rate.lims) | max(rate.range) < max(rate.lims)) {
			warning(paste("Did not override rate.lims: the specified rate.range (", rate.range[1], ", ", rate.range[2], ") did not contain all the values for the observed rates (", rate.lims[1], ", ", rate.lims[2], ")"))
		} else {
			rate.lims <- rate.range
		}
	}
	rate.lims[1] <- rate.lims[1] - lims.percentage.correction*abs(rate.lims[1])
	rate.lims[2] <- rate.lims[2] + lims.percentage.correction*abs(rate.lims[2])
	
	rate.tree <- contMapGivenAnc(tree= tree.to.plot, x=rates.tips, plot=FALSE, anc.states=rates.internal, lims=rate.lims, ...)
	#change colors
	rate.colors <- colorRampPalette(rate.colors, space="Lab")(length(rate.tree$cols))
	rate.tree$cols[] <- rate.colors
	#rate.tree$cols[] <- adjustcolor(rate.tree$cols[], alpha.f=0.3)
	plot(rate.tree, outline=FALSE, lwd=edge.width.rate, legend=FALSE, type=type, fsize=fsize, ...)
	par(fig=c(0,1, 0, 1), new=TRUE)
	#state.tree <- contMapGivenAnc(tree=hisse.object$phy, x=ConvertToBinaryState(hisse.object$tip.mat, state.0.indices=state.0.indices), plot=FALSE, anc.states=ConvertToBinaryState(hisse.object$node.mat, state.0.indices=state.0.indices))
	
	state.lims <- range(c(states.tips, states.internal))
	state.lims[1] <- state.lims[1] - lims.percentage.correction*abs(state.lims[1])
	state.lims[2] <- state.lims[2] + lims.percentage.correction*abs(state.lims[2])

	state.tree <- contMapGivenAnc(tree=tree.to.plot, x=states.tips, plot=FALSE, anc.states=states.internal, lims=state.lims, ...)
	#state.colors <- grey(seq(1,0,length.out=length(state.tree$cols)))
	state.colors <- colorRampPalette(state.colors, space="Lab")(length(rate.tree$cols))
	state.tree$cols[]<- state.colors
	plot(state.tree, outline=FALSE, lwd=edge.width.state, legend=FALSE, type=type, fsize=fsize, ...)
	
	if(legend!="none") {
		par(fig=legend.position, new=TRUE)
		plot(x=c(-0.1,1.1), y=c(-1.5,1.5), xlab="", ylab="", bty="n", type="n", xaxt="n", yaxt="n")
		rect(-0.1,-1.1,1.1,1.1, border=NA, col=legend.bg)
		rates.to.plot <- c()
		states.to.plot <- c()
		if(legend=="all" | legend=="tips") {
			rates.to.plot <- append(rates.to.plot, rates.tips)
			states.to.plot <- append(states.to.plot, states.tips)
		}
		if(legend=="all" | legend=="internal") {
			rates.to.plot <- append(rates.to.plot, rates.internal)
			states.to.plot <- append(states.to.plot, states.internal)
		}
		
		if(legend.kernel.rates=="auto") {
			if(length(unique(rates.to.plot))<=4) {
				legend.kernel.rates <- "hist"	
			} else {
				legend.kernel.rates <- "rectangular"	
			}
		}
		if(legend.kernel.states=="auto") {
			if(length(unique(states.to.plot))<=4) {
				legend.kernel.states <- "hist"	
			} else {
				legend.kernel.states <- "rectangular"	
			}
		}

		rates.density <- GetNormalizedDensityPlot(rates.to.plot, rate.lims, legend.kernel.rates)
		states.density <- GetNormalizedDensityPlot(states.to.plot, state.lims, legend.kernel.states)
		states.density$y <- (-1) * states.density$y #so it gets drawn below the other one
		# rates.density <- c()
		# states.density <- c()
		# if (legend.kernel=="hist") {
			# rates.density<-hist(rates.to.plot, breaks=seq(from=rate.lims[1], to=rate.lims[2], length.out = max(100,nclass.Sturges(rates.to.plot)+2)), plot=FALSE)
			# states.density<-hist(states.to.plot, breaks=seq(from=state.lims[1], to=state.lims[2], length.out = max(100,nclass.Sturges(states.to.plot)+2)), plot=FALSE)
			# rates.density$x <- rates.density$mid
			# rates.density$y <- rates.density$density
			# states.density$x <- states.density$mid
			# states.density$y <- states.density$density
		# } else {
			# rates.density <- density(rates.to.plot, from=rate.lims[1], to=rate.lims[2], kernel=legend.kernel)
			# states.density <- density(states.to.plot, from=state.lims[1], to=state.lims[2], kernel=legend.kernel)
		# }
		# rates.density$x <- (rates.density$x - rate.lims[1]) / (rate.lims[2]-rate.lims[1]) #so it goes from zero to one
		# rates.density$y <- rates.density$y/max(rates.density$y)
		# states.density$y <- (-1) * states.density$y/max(states.density$y)
		# states.density$x <- (states.density$x - state.lims[1]) / (state.lims[2]-state.lims[1]) #so it goes from zero to one

		# if(legend=="rect") {
			# rates.density$y <- rep(1, length(rates.density$y))
			# states.density$y <- rep(-1, length(states.density$y))
		# }
		par(lend=1)
		segments(x0=rates.density$x, y0=rep(0, length(rates.density$y)), y1=rates.density$y, col=rate.colors[1+as.integer(round((length(rate.colors)-1)* rates.density$x))], lwd=ifelse(legend.kernel.rates=="hist",4,1))
		text(x=0, y=1.2, labels=format(rate.lims[1], digits=2), cex=legend.cex)
		text(x=1, y=1.2, labels=format(rate.lims[2], digits=2), cex=legend.cex) 
		text(x=0.5, y=1.2, labels=rate.param, cex=legend.cex)	
#		lines(rates.density$x, rates.density$y, lwd=0.5, col="gray")
		
		segments(x0=states.density$x, y0=rep(0, length(states.density$y)), y1=states.density$y, col=state.colors[1+as.integer(round((length(state.colors)-1)* states.density$x))], lwd=ifelse(legend.kernel.states=="hist",4,1))
		text(x=0, y=-1.2, labels="0", cex=legend.cex)
		text(x=1, y=-1.2, labels="1", cex=legend.cex) 
		text(x=0.5, y=-1.2, labels="State", cex=legend.cex)	
#		lines(states.density$x, states.density$y, lwd=0.5, col="gray")
#		lines(rates.density$x, 0*rates.density$y, lwd=0.5, col="gray")
	
	}
	
	return(list(rate.tree=rate.tree, state.tree=state.tree))
}

GetNormalizedDensityPlot <- function(x, limits, kernel, min.breaks=100) {
	x.density <- c()
	
	if (kernel=="hist") {
		x.density<-hist(x, breaks=seq(from= limits[1], to= limits[2], length.out = max(min.breaks,nclass.Sturges(x)+2)), plot=FALSE)
		x.density$x <- x.density$mid
		x.density$y <- x.density$density
		x.density$x <- x.density$x[which(x.density$y>0)] #since the line is thick, do not plot it if zero
		x.density$y <- x.density$y[which(x.density$y>0)]
	} else {
		if(kernel=="traditional") {
			x.density <- density(x, from=limits[1], to=limits[2])
			x.density$y <- rep(1, length(x.density$y))
		} else {
			x.density <- density(x, from=limits[1], to=limits[2], kernel=kernel)
		}
	}	
	x.density$x <- (x.density$x - limits[1]) / (limits[2]-limits[1]) #so it goes from zero to one
	if(kernel!="traditional") {
		x.density$y <- x.density$y/max(x.density$y)
	}
	return(x.density)
}

# The following function is from phytools, by Liam Revell, which is released under GPL2+
# Modified by Brian O''Meara, June 9, 2015

# function plots reconstructed values for ancestral characters along the edges of the tree
# written by Liam J. Revell 2012, 2013, 2014, 2015

contMapGivenAnc <-function(tree,x,res=100,fsize=NULL,ftype=NULL,lwd=4,legend=NULL,
lims=NULL,outline=TRUE,sig=3,type="phylogram",direction="rightwards",
plot=TRUE,anc.states=NULL,...){
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(method)) method<-list(...)$method
	else method<-"fastAnc"
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	h<-max(nodeHeights(tree))
	steps<-0:res/res*max(h)
	H<-nodeHeights(tree)
	a <- anc.states #BCO modified
	if(is.null(a)) { #BCO put this in if loop
		if(method=="fastAnc") a<-fastAnc(tree,x) 
		else { 
			fit<-anc.ML(tree,x)
			a<-fit$ace
			if(!is.null(fit$missing.x)) x<-c(x,fit$missing.x)
		}
	} #end BCO if loop
	names(x) <- tree$tip.label[as.numeric(names(x))]
	y<-c(a,x[tree$tip.label]); names(y)[1:length(tree$tip)+tree$Nnode]<-1:length(tree$tip)
	A<-matrix(y[as.character(tree$edge)],nrow(tree$edge),ncol(tree$edge))
	cols<-rainbow(1001,start=0,end=0.7); names(cols)<-0:1000
	if(is.null(lims)) lims<-c(min(c(a,x)),max(c(a,x))) #modified by BCO to include anc state in range for lims
	trans<-0:1000/1000*(lims[2]-lims[1])+lims[1]; names(trans)<-0:1000
	tree$maps <- list(rep(rep(NA, 2), nrow(tree$edge)))
	for(i in 1:nrow(tree$edge)){
		XX<-cbind(c(H[i,1],steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))]),
				  c(steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))],H[i,2]))-H[i,1]
		YY<-rowMeans(XX)
		if(!all(YY==0)){
			b<-vector()
			for(j in 1:length(YY))
			b[j]<-(A[i,1]/YY[j]+A[i,2]/(max(XX)-YY[j]))/(1/YY[j]+1/(max(XX)-YY[j]))
		} else b<-A[i,1]
		d<-sapply(b,getState,trans=trans)
		tree$maps[[i]]<-XX[,2]-XX[,1]
		names(tree$maps[[i]])<-d
	}
	tree$mapped.edge<-makeMappedEdge(tree$edge,tree$maps)
	tree$mapped.edge<-tree$mapped.edge[,order(as.numeric(colnames(tree$mapped.edge)))]
	xx<-list(tree=tree,cols=cols,lims=lims)
	class(xx)<-"contMap"
	if(plot) plot.contMap(xx,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,
						  sig=sig,type=type,mar=mar,direction=direction,offset=offset,hold=hold)
	invisible(xx)
}


# The following function is from phytools, by Liam Revell, which is released under GPL2+
# It is not exported from phytools, so the options are phytools:::getState or copy it here
# Given that phytools could change non-exported functions in a way that breaks the above code
# I have elected to copy it unchanged.
# Brian O'Meara, June 10, 2015

# function
# written by Liam J. Revell 2012
getState<-function(x,trans){
	i<-1
	state <- names(trans)[1] #BCO: added to prevent error when while loop not entered
	while(x>trans[i]){
		state<-names(trans)[i]
		i<-i+1
	}
	return(state)
}

# The following function is from phytools, by Liam Revell, which is released under GPL2+
# It is not exported from phytools, so the options are phytools:::makeMappedEdge or copy it here
# Given that phytools could change non-exported functions in a way that breaks the above code
# I have elected to copy it unchanged. Phytools remains a 
# Brian O'Meara, June 10, 2015

# make a mapped edge matrix
# written by Liam J. Revell 2013

makeMappedEdge<-function(edge,maps){
	st<-sort(unique(unlist(sapply(maps,function(x) names(x)))))
	mapped.edge<-matrix(0,nrow(edge),length(st))
	rownames(mapped.edge)<-apply(edge,1,function(x) paste(x,collapse=","))
	colnames(mapped.edge)<-st
	for(i in 1:length(maps)) 
		for(j in 1:length(maps[[i]])) 
			mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
	return(mapped.edge)
}

ConvertToBinaryState <- function(x) {
	x.trimmed <- x[,-1]
	state.0.indices <- which(grepl("0", colnames(x.trimmed)))
	state.1.indices <- which(grepl("1", colnames(x.trimmed)))
	x0.trimmed <- x.trimmed[, state.0.indices]
	if (is.null(dim(x0.trimmed)[1])) { #is a vector
		x0.trimmed <- matrix(data=x0.trimmed, nrow=length(x0.trimmed), ncol=1)
	}
	x1.trimmed <- x.trimmed[, state.1.indices]
	if (is.null(dim(x1.trimmed)[1])) { #is a vector
		x1.trimmed <- matrix(data=x1.trimmed, nrow=length(x1.trimmed), ncol=1)
	}
	result.vector.0 <- apply(x0.trimmed, 1, sum)
	result.vector.1 <- apply(x1.trimmed, 1, sum)
	result.vector <- result.vector.0 / (result.vector.0 + result.vector.1)
	names(result.vector) <- x[,1]
	return(result.vector)
}

ConvertToRate <- function(x, rate.vector) {
	x.trimmed <- x[,-1]
	x.trimmed <- x.trimmed / rowSums(x.trimmed) #normalize to 1 (which it should be already)
	result.vector <- rep(0, dim(x.trimmed)[1])
	for (i in sequence(dim(x.trimmed)[2])) {
		result.vector <- result.vector + rate.vector[i] * x.trimmed[,i]	#get weighted mean
	}
	names(result.vector) <- x[,1]
	return(result.vector)
}

ConvertManyToRate <- function(hisse.results, rate.param, which.element) {
	AIC.weights <- GetAICWeights(hisse.results)
	storage.matrix <- matrix(nrow=dim(hisse.results[[1]][[which.element]])[1], ncol=0)
	for (i in sequence(length(hisse.results))) {
		rate.vector <- hisse.results[[i]]$rates.mat[rate.param,]
		storage.matrix <- cbind(storage.matrix, ConvertToRate(x=hisse.results[[i]][[which.element]], rate.vector=hisse.results[[i]]$rates.mat[rate.param,]))	
	}
	final.results <- apply(storage.matrix, 1, weighted.mean, w=AIC.weights)
	return(final.results)
}

ConvertManyToBinaryState <- function(hisse.results, which.element) {
	AIC.weights <- GetAICWeights(hisse.results)
	storage.matrix <- matrix(nrow=dim(hisse.results[[1]][[which.element]])[1], ncol=0)
	for (i in sequence(length(hisse.results))) {
		storage.matrix <- cbind(storage.matrix, ConvertToBinaryState(x=hisse.results[[i]][[which.element]]))	
	}
	final.results <- apply(storage.matrix, 1, weighted.mean, w=AIC.weights)
	return(final.results)
}

GetAICWeights <- function(hisse.results) {
	AIC.vector <- sapply(hisse.results, "[[", "aic")
	delta.AIC.vector <- AIC.vector - min(AIC.vector)
	rel.likelihood <- exp(-0.5 * delta.AIC.vector)
	AIC.weight.vector <- rel.likelihood / sum(rel.likelihood)
	return(AIC.weight.vector)
}


GetRateRange <- function(x, rate.param) {
	hisse.results <- x
	if(class(hisse.results)=="hisse.states") { #we have to make a list so we can run this generally
		tmp.list <- list()
		tmp.list[[1]] <- hisse.results
		hisse.results <- tmp.list	
	}
	all.rates <- sapply(hisse.results, "[[", "rates.mat", simplify=FALSE)
	return(range(unname(unlist(sapply(all.rates, GetRelevantRowEntries, rate.param="turnover")))))	
}

GetRelevantRowEntries <- function(x, rate.param) {
	return(x[which(rownames(x)==rate.param),])
}
