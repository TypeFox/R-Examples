# utils.R
# 
# Other useful methods for HMMs
#
# Author: Peter Humburg
###############################################################################

############################ auxiliary functions ##############################
## functions to handle the log space computations
## calculate the log of the sum from the log of the summands
## base: base of logarithm, 0 indicates natural log
logSum <- function(x,y=NULL,base=0){
	if(is.null(y)){
		if(is.vector(x)){
			lsum <- x[1]
			if(length(x) == 1) return(lsum)
			for(i in 2:length(x)){
				lsum <- logSum(lsum,x[i]) 
			}
			return(lsum)
		}
		if(is.matrix(x)){
			lsum <- apply(x,1,logSum)
			lsum <- logSum(lsum)
			return(lsum)
		}
	}
	if(is.matrix(x) && is.matrix(y)){
		if(sum(dim(x) != dim(y))){
			stop("logSum: Matrix dimensions do not agree.")
		}
		lsum <- x
		for(i in 1:dim(x)[1]){
			for(j in 1:dim(x)[2]){
				lsum[i,j] <- logSum(x[i,j],y[i,j])
			}
		}
		return(lsum)
	}
	if(length(x) > 1 || length(y) > 1){
		if(length(x) != length(y)){
			stop("logSum: vectors of different length.")
		}
		lsum <- numeric(length(x))
		for(i in 1:length(x)){
			lsum[i] <- logSum(x[i],y[i])
		}
		return(lsum)
	}
	
	.Call("_log_sum",x,y,as.integer(base))
}

## Use combination of Viterbi training and Baum-Welch algorithm to fit HMM
viterbiEM <- function(hmm,data,max.iter=c(5,15),eps=0.01,verbose=0,...){
	args <- list(...)
	if(is.null(args$viterbi)) vit.args <- args[names(args) != "df"] else vit.args <- args$viterbi
	if(is.null(args$baumWelch)) bw.args <- args[names(args) != "df"] else bw.args <- args$baumWelch
	df <- args$df
	
	max.iter <- rep(max.iter, length = 2)
	eps <- rep(eps, length = 2)
	
	if(class(df) == "list"){
		df.vit <- df[[1]]
		df.bw <- df[[2]]
	} else{
		df.vit <- df
		df.bw <- df
	}
	if(verbose >=1) message("Viterbi training: ", max.iter[1], " iterations\n")
	hmm.vit <- do.call("viterbiTraining",c(list(hmm,data,max.iter=max.iter[1],eps=eps[1],
	             df=df.vit,verbose=verbose-1),vit.args))
	
	## check that all states are accessible
	## get prior
	if(is.null(bw.args$trans.prior)){
		## no prior
		bw.args$trans.prior <- matrix(0,ncol=length(hmm),nrow=length(hmm))
	} else if(is.logical(bw.args$trans.prior) && bw.args$trans.prior){
		bw.args$trans.prior <- hmm@transition.matrix
	}
	if(!is.matrix(bw.args$trans.prior) || dim(bw.args$trans.prior) != dim(hmm@transition.matrix)){
		stop("Illegal prior transition distribution")
	}
	## check prior initial state distribution
	if(is.null(bw.args$init.prior)){
		## no prior
		bw.args$init.prior <- numeric(length(hmm))
	} else if(is.logical(bw.args$init.prior) && bw.args$init.prior){
		bw.args$init.prior <- hmm@init
	} 
	if(!is.vector(bw.args$init.prior) || length(bw.args$init.prior) != length(hmm)){
		stop("Illegal prior initial state distribution")
	}
	access.trans <- hmm.vit@transition.matrix + bw.args$trans.prior
	access.init <- hmm.vit@init + bw.args$init.prior
	diag(access.trans) <- access.init
	access <- apply(access.trans, 2, sum) > 0
	if(sum(access == 0)){
		warning("Detected ", sum(access == 0), " state(s) without inbound transitions (",
				paste(states(hmm.vit)[which(access == 0)], collapse=', '), ").\n",
				"Continuing with uniform prior.", call.=FALSE, immediate.=TRUE)
		bw.args$trans.prior <- bw.args$trans.prior + 1/length(hmm.vit)
		bw.args$init.prior <- bw.args$init.prior + 1/length(hmm.vit)
		## need to update current transition probabilities
		init <- .baumWelchInit(hmm.vit, data)
		trans <- .baumWelchTransition(hmm.vit, data, init$alpha, init$beta, bw.args$trans.prior, bw.args$init.prior)
		hmm.vit@transition.matrix <- trans$transition
	}
	
	if(verbose >=1) message("EM algorithm: ", max.iter[2], " iterations\n")
	hmm.bw <- do.call("baumWelch",c(list(hmm.vit,data,max.iter=max.iter[2],
				eps=eps[2],df=df.bw, verbose=verbose-1),bw.args))
				
	hmm.bw
}

## calculate posterior probability for each state of hmm
posterior <- function(data, hmm, log=TRUE){
	alpha <- forward(hmm, data)$alpha.scaled
	beta <- backward(hmm, data)
	
	post <- alpha + beta
	post <- t(t(post) - apply(post,2,logSum))
	if(!log) post <- exp(post)
	return(post)
}

## takes a matrix containing information about enriched regions, a data frame with information about probe positions,
## and a vector with posterior probabilities for each probe
## returns a data frame in gff format invisibly.
## output is saved to 'file' if a file name or connection is provided
reg2gff <- function(regions,post,probe.pos,file=NULL,score.fun=mean,source="tHMM",feature.type="posterior_prob",class="ChIP_region",name="tHMM"){
	chr <- probe.pos[regions[1,],"chromosome"]
	score <- apply(regions,2,function(reg,p)score.fun(p[reg[1]:reg[2]]),post)
	start <- probe.pos[regions[1,],"position"]
	end <- probe.pos[regions[2,],"position"]
	gff <- data.frame(chr=chr,source=I(source),type=I(feature.type),start=start,end=end,
						score=score,strand=I('.'),phase=I('.'),group=I(paste(class,paste(name,1:(dim(regions)[2]),sep='_'))))
	if(!is.null(file)){
		write.table(gff,file=file,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
	}
	
	return(invisible(gff))
}

## takes data from a gff file and probe positions, returns a logical vector indicating probes in
## annotated regions
## If gff is a character string it is assumed to be the name of a gff file which is read and used subsequently
gff2index <- function(gff,pos){
	if(is.character(gff) || is(gff,"connection")){
		gff <- read.delim(gff,comment.char='#',header=FALSE)
	}
	index <- logical(dim(pos)[1])
	for(i in 1:dim(gff)[1]){
		chr <- gff[i,1]
		start <- gff[i,4]
		end <- gff[i,5]
		index <- index | (pos[,1] == chr & pos[,2] >= start & pos[,2] <= end)
	}
	index
}

## takes a logical vector indicating positive and negative probes
## returns a list with components 'positive' and 'negative' providing length
## information for positive and negative regions
region.length <- function(probes, min.len=1){
	positive <- numeric()
	negative <- numeric()
	pos.count <- as.numeric(probes[1])
	neg.count <- 0
	if(!probes[1]) neg.count <- 1
	
	for(p in probes[-1]){
		## extend current region
		if(p & pos.count > 0) pos.count <- pos.count + 1
		if(!p & neg.count > 0) neg.count <- neg.count + 1
		## start of new region
		if(p & pos.count == 0){
			if(neg.count >= min.len) negative <- c(negative,neg.count)
			neg.count <- 0
			pos.count <- 1
		}
		if(!p & neg.count == 0){
			if(pos.count >= min.len) positive <- c(positive,pos.count)
			pos.count <- 0
			neg.count <- 1
		}
	}
	## add last region
	if(pos.count > min.len) positive <- c(positive,pos.count)
	if(neg.count > min.len) negative <- c(negative,neg.count)
	
	ret <- list()
	ret[["positive"]] <- positive
	ret[["negative"]] <- negative
	
	ret
}

## Takes a vector containing calls for each probe, returns a matrix with
## start and end points of the regions of interest in its first and second row respectively. 
region.position <- function(probe.calls,region=TRUE){
	probe.index <- probe.calls == region
	idx <- 1:length(probe.calls)
	probe.index2 <- idx[probe.index]
	probe.dist <- diff(probe.index2)
	region.idx <- probe.dist > 1
	idx <- 1:length(probe.dist)
	region.idx2 <- idx[region.idx]
	end <- c(probe.index2[region.idx2],probe.index2[length(probe.index2)])
	start <- c(probe.index2[1],probe.index2[region.idx2+1])
	
	rbind(start,end)
}

## takes a matrix with information about enriched regions (in terms of probe indices),
## posterior probabilities for each probe and a data frame with probe positions
## removes regions that are shorter than min.length and have a score of less than min.score
## the score for each region is calculated by applying summary.fun to the posterior probabilities 
## of all probes in the region
remove.short <- function(regions, post, probe.pos, min.length=1000, min.score=0.8, summary.fun=mean){
	regions.pos <- matrix(probe.pos[regions, "position"], nrow = 2, ncol = dim(regions)[2])
	reg.len <- apply(regions.pos, 2, diff)
	if(max(post) < 0) post <- exp(post)
	reg.score <- apply(regions, 2, function(reg) summary.fun(post[reg[1]:reg[2]]))
	reg.idx <- reg.len >= min.length | reg.score >= min.score
	regions[ , reg.idx]
}

## generate simulated data set by sampling from real data using results of TileMap analysis
## data is a data.frame containing probe positions and measurements
## group is either a logical vector indicating the position of positive probes or
##     the name of a gff file containing information about positive regions
## the number of positive regions per sequence is between min(pos.range) and max(pos.range)
##     for each sequence the actual number is sampled from a uniform distribution
## num.seq sequences are generated
## probe positions are generated using a distance of gap between probes of the same sequence
##      and split.gap between different sequences
## returns a data matrix
generate.data <- function(data,group,pos.range=c(1,10),num.seq=100,gap=35,split.gap=1000,min.len=2){
	## read region information from file if necessary
	if(is.character(group) || is(group,"connection")){
		gff <- read.delim(group,comment.char='#',header=FALSE)
		group <- gff2index(gff,data[,c(1,2)])
	}
	## get length distributions
	reg.len <- region.length(group,min.len=min.len)
	
	## generate observation sequences
	pos.count <- runif(num.seq,min=min(pos.range),max=max(pos.range))
	probe.pos <- 1
	output <- matrix(ncol=dim(data)[2],nrow=0)
	index <- 1:dim(data)[1]
	pos.prob <- sum(group)/length(group)
	start <- sample(c(0,1),num.seq,prob=c(1-pos.prob,pos.prob),replace=TRUE)
	end <- sample(c(0,1),num.seq,prob=c(1-pos.prob,pos.prob),replace=TRUE)
	state.seq <- list()
	for(i in 1:num.seq){
		sequence <- matrix(ncol=dim(data)[2]-2,nrow=0)
		state.seq[[i]] <- logical()
		if(pos.count[i] > 0){
			if(start[i]){
				## positive region
				seq.len <- sample(reg.len$positive,1)
				seq.idx <- sample(index[group],seq.len)
				sequence <- rbind(sequence,data[seq.idx,c(-1,-2)])
				pos.count[i] <- pos.count[i]-1
				state.seq[[i]] <- c(state.seq[[i]],rep(TRUE,each=seq.len))
			}
			while(pos.count[i] > end[i]){
				## negative region
				seq.len <- sample(reg.len$negative,1)
				seq.idx <- sample(index[!group],seq.len)
				sequence <- rbind(sequence,data[seq.idx,c(-1,-2)])
				state.seq[[i]] <- c(state.seq[[i]],rep(FALSE,each=seq.len))
				## positive region
				seq.len <- sample(reg.len$positive,1)
				seq.idx <- sample(index[group],seq.len)
				sequence <- rbind(sequence,data[seq.idx,c(-1,-2)])
				pos.count[i] <- pos.count[i]-1
				state.seq[[i]] <- c(state.seq[[i]],rep(TRUE,each=seq.len))
			}
			## negative region
			seq.len <- sample(reg.len$negative,1)
			seq.idx <- sample(index[!group],seq.len)
			sequence <- rbind(sequence,data[seq.idx,c(-1,-2)])
			state.seq[[i]] <- c(state.seq[[i]],rep(FALSE,each=seq.len))
			if(end[i]){
				## positive region
				seq.len <- sample(reg.len$positive,1)
				seq.idx <- sample(index[group],seq.len)
				sequence <- rbind(sequence,data[seq.idx,c(-1,-2)])
				state.seq[[i]] <- c(state.seq[[i]],rep(TRUE,each=seq.len))
			}
		} else{
			## negative region
			seq.len <- sample(reg.len$negative,1)
			seq.idx <- sample(index[!group],seq.len)
			sequence <- rbind(sequence,data[seq.idx,c(-1,-2)])
			state.seq[[i]] <- c(state.seq[[i]],rep(FALSE,each=seq.len))
		}
		## generate probe positions
		positions <- seq(from=probe.pos,by=gap,length.out=dim(sequence)[1])
		chr <- rep("chr1",each=dim(sequence)[1])
		output <- rbind(output,cbind(chr,positions,sequence))
		probe.pos <- positions[length(positions)]+split.gap
	}
	ret <- list()
	ret[["observation"]] <- output
	ret[["regions"]] <- state.seq
	ret
} 

## get initial partition of data into states
.get.init <- function(data,nstate=2,iter.max=20,nstart=3){
	cluster <- kmeans(data,nstate,iter.max=iter.max,nstart=nstart)
	means <- cluster$centers
	vars <- cluster$withinss/(cluster$size - 1)
	
	ret <- list()
	ret[["mean"]] <- means
	ret[["var"]] <- vars
	ret[["cluster"]] <- cluster$cluster
	ret[["size"]] <- cluster$size
	ret	
}

## setting up an HMM with initial estimates from the data
hmm.setup <- function(data,state=c("enriched","non-enriched"),probe.region=35,frag.size=1000,pos.state=1,
						em.type="tDist",max.prob=1,df=9){
	comp <- numeric(length(state)) + 1
	names(comp) <- state
	## cluster data to obtain initial estimates
	if(is(data,"list")) data <- c(data,recursive=TRUE)
	cl <- .get.init(data,nstate=c(min(data),max(data)))
	
	mean.pos <- cl$mean[1]
	mean.neg <- cl$mean[2]
	var.pos <- cl$var[1]
	var.neg <- cl$var[2]
	weight.pos <- 1
	weight.neg <- 1
	emission.comp <- list()
	emission <- list()

	emission.comp[[state[1]]] <- cbind(weight.pos,mean.pos,var.pos)
	emission.comp[[state[2]]] <- cbind(weight.neg,mean.neg,var.neg)
	for(i in 1:length(state)){
		if(em.type == "tDist"){
			emission[[state[i]]] <- new("tDist",mean=emission.comp[[state[i]]][,2],
									var=emission.comp[[state[i]]][,3],df=df[(i%%length(df))+1])
		} else{
			stop(paste("Requested emission distribution of type",em.type,"is not supported."))	
		}
	}

	## set initial transition probabilities according to average fragment size
	## we don't know how many enriched regions to expect, but the number 
	## of sufficiently extreme observations should be a good starting point
	pos.index <- cl$cluster == 1
	neg.index <- !pos.index
	freq.pos <- sum(pos.index)/length(data)
	freq.pos <- freq.pos*(probe.region/frag.size)
	
	if(pos.state == 1){
		prob.aa <- min(1-probe.region/frag.size,max.prob)
		prob.bb <- min(1-freq.pos,max.prob)
	} else{
		prob.aa <- min(1-freq.pos,max.prob)
		prob.bb <- min(1-probe.region/frag.size,max.prob)
	}
	
	neg.state <- 2
	if(pos.state == 2) neg.state <- 1
	transition <- list()
	transition[[state[1]]] <- new("discDist",alpha=state,prob=c(prob.aa,1-prob.aa))
	transition[[state[2]]] <- new("discDist",alpha=state,prob=c(1-prob.bb,prob.bb))
	
	## choose initial distribution as stationary distribution of the MC with above transition matrix
	a <- transition[[state[1]]][2]
	b <- transition[[state[2]]][1]
	pi <- new("discDist",alpha=state,prob=c(b/(a+b),a/(a+b)))
	
	## create HMM with initial parameters
	hmm.init <- new("contHMM",transition=transition,emission=emission,init=pi)
	hmm.init
}

## a more convenient way to create an HMM
getHMM <- function(params,snames){
	trans <- list()
	trans[[snames[1]]] <- new("discDist",alpha=snames,prob=c(1-params$a[1],params$a[1]))
	trans[[snames[2]]] <- new("discDist",alpha=snames,prob=c(params$a[2],1-params$a[2]))

	emiss <- list()
	emiss[[snames[1]]] <- new("tDist",mean=params$mu[1],var=params$sigma[1],df=params$nu[1])
	emiss[[snames[2]]] <- new("tDist",mean=params$mu[2],var=params$sigma[2],df=params$nu[2])

	init <- new("discDist",alpha=snames,prob=c(params$a[2]/(params$a[2]+params$a[1]),
		params$a[1]/(params$a[2]+params$a[1])))

	new("contHMM",transition=trans,emission=emiss,init=init)
}

## calculate 'shrinkage t' statistic
shrinkt.st <- function(X,L,h0.mean=0,...){
	## ensure that there are at least two samples per group
	tbl <- table(L)
	if(min(tbl) == 1){
		stop("Group ", names(tbl)[which.min(tbl)], " contains only one sample. At least two samples are required per group.")
	}
	
	if(max(L) == 1){	
		shrink.var <- var.shrink(t(X))
		means <- apply(X,1,mean)
		t.stat <- mapply(function(m,v,h0,n) (m - h0)/sqrt(v/n),means,shrink.var,
			MoreArgs=list(h0.mean,length(L),...))
	} else if(max(L) == 2){
		t.stat <- shrinkt.stat(t(X),L,...)
	} else stop(paste("Illegal number of groups! Expected 1 or 2, found", max(L)))
	
	t.stat	
}