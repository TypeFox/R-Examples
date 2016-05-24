tsvq <-
function(x, K=nrow(x), row.labs = 1:nrow(x),ntry=20,verbose=FALSE,
as.hclust=TRUE,trace=FALSE)
{
	if (is.data.frame(x)) x <- as.matrix(x)
	if (is.null(row.names(x))) row.names(x) <- 1:nrow(x)
	trytwomeans <- function(x,centers,iter.max=10,ntry){
		if (nrow(centers)!=2) stop(paste('trytwomeans: trying to split into',
			nrow(centers), 'groups instead of 2\n'))
		best <- kmeans(x,centers,iter.max)
		if (ntry>1){
			for (i in 2:ntry){
				seed <- sample(nrow(x),nrow(centers),replace=FALSE)
      	while (all(x[seed[1],]==x[seed[2],],na.rm=TRUE))
					seed <- sample(nrow(x),nrow(centers),replace=FALSE)
				temp <- kmeans(x,x[seed,],iter.max)
				if (sum(temp$withinss) < sum(best$withinss)){
					best <- temp
				}
			}
		}
		best
	}

	center <- apply(x,2,mean,na.rm=TRUE)
	interior.list <- list(list(rlabs=row.labs,center=center,size=nrow(x),
		wss=sum(scale(x,scale=FALSE)^2),path=''))
	if (verbose)
		cat("Start WSS:", format(round(interior.list[[1]]$wss, 2)), "\n")
	means2 <- function(x, rlabs,verbose)
	{
		if(length(rlabs) == 2) {
			list(wss = c(0, 0), pair = list(list(rlabs = rlabs[1], 
				center = x[1,  ], size = 1, wss = 0), list(
				rlabs = rlabs[2], center = x[2,  ], size = 1, 
				wss = 0)))
		}
		else {
			if(!is.matrix(x) || is.null(x) || nrow(x)==0) browser()
			jseq <- seq(nrow(x))
			seed <- sample(jseq, 2, replace = FALSE)
			# next loop makes sure we don't have two duplicate points
			while (all(x[seed[1],]==x[seed[2],],na.rm=TRUE))
				seed <- sample(jseq, 2, replace = FALSE)
			#cat('\n',length(interior.list),'seed=',seed,':  ')
			fit0 <- trytwomeans(x, x[seed,  ],ntry=ntry)
			junk <- as.list(1:2)
			wss <- fit0$withinss
			if(verbose) cat("Split pair WSS:", format(round(wss, 2)), "\n")
			for(i in 1:2) {
				junk[[i]] <- list(rlabs = rlabs[fit0$cluster == 
				  i], center = fit0$centers[i,  ], size = fit0$
				  size[i], wss = wss[i])
			}
			o <- order(wss)
			list(wss = wss[o], pair = junk[o])
		}
	}
	junk <- means2(x, row.labs,verbose)
	wss <- junk$wss
	tree.list <- junk$pair
	tree.list[[1]]$path <- "Left"
	tree.list[[2]]$path <- "Right"
	k <- 2
	while(k < K) {
		if (!verbose & trace) cat(k)
		if (verbose) cat("Split ", k, "WSS before", format(round(wss[k], 2)), "\n")
		#else cat(k)  # uncomment to get some indication of how tsvq is progressing
		rlabs <- tree.list[[k]]$rlabs
		
		#cat('starting means2...')
		junk <- means2(x[rlabs,  ], rlabs,verbose)
		#cat('...and done\n')
		path <- tree.list[[k]]$path
		jpair <- junk$pair
		jpair[[1]]$path <- c(path, "Left")
		jpair[[2]]$path <- c(path, "Right")
		interior.list[[length(interior.list)+1]] <- tree.list[k][[1]]
		wss <- c(wss[ - k], junk$wss)
		tree.list <- c(tree.list[ - k], jpair)
		o <- order(wss)
		wss <- wss[o]
		tree.list <- tree.list[o]
		k <- k + 1
	}
	result <- list(tree.list=tree.list,interior.list=interior.list)
	class(result) <- 'tsvq'
	if (as.hclust) return(tsvq2hclust(result))
	if (verbose) cat('\n')
	else return(result)
}

redo.fused.tsvq <- function(obj,x,...)
{
	newobj <- obj
	# recalculate the wss for the interior nodes of the tree
	for (i in 1:length(obj[[2]]))
	{
		if (length(obj[[2]][[i]]$rlabs)>1)
		{
			xx <- x[obj[[2]][[i]]$rlabs,]
			center <- apply(xx,2,mean,na.rm=TRUE)
			newobj[[2]][[i]]$wss <- sum(scale(xx,scale=FALSE)^2)
		}
	}
	# for each fused group, redo tsvq and update list.
	for (i in 1:length(obj[[1]]))
	{
		if (length(obj[[1]][[i]]$rlabs)>1)
		{
			rl <- obj[[1]][[i]]$rlabs
			xx <- x[rl,]
			junk <- tsvq(xx,length(rl),row.labs=1:length(rl),as.hclust=FALSE,...)
			newobj[[1]][[i]] <- 'kill'
			for (j in 1:length(junk[[1]]))
				junk[[1]][[j]]$rlabs <- rl[junk[[1]][[j]]$rlabs]
			for (j in 1:length(junk[[2]]))
				junk[[2]][[j]]$rlabs <- rl[junk[[2]][[j]]$rlabs]
			newobj[[1]][1:length(junk[[1]])+length(newobj[[1]])] <- junk[[1]]
			newobj[[2]][1:length(junk[[2]])+length(newobj[[2]])] <- junk[[2]]
		}
	}
	i <- 1
	while (i < length(newobj[[1]]))
	{
		if (is.character(newobj[[1]][[i]])) newobj[[1]][[i]] <- NULL
		else i <- i+1
	}
	newobj
}

tsvq2hclust <- function(obj){
	obj <- obj$interior.list
	descendants <- function(m,k){
  # the done object indicates what rows of m were used
    done <- k
    if (m[k,1]<0) left <- -m[k,1]
    else {
      junk <- descendants(m,m[k,1])
      left <- junk[[1]]
      done <- c(done,junk[[2]])
    }
    if (m[k,2]<0) right <- -m[k,2]
    else {
      junk <- descendants(m,m[k,2])
      right <- junk[[1]]
      done <- c(done,junk[[2]])
    }
    return(list(c(left,right),done))
  }
	n <- length(obj)+1
	height <- rep(0,n-1)
	merge <- matrix(0,n-1,2)
	for (i in 1:(n-1)) height[i] <- obj[[i]]$wss
	hord <- order(height)
	height <- height[hord]
	obslist <- vector('list',n-1)
	for (i in 1:(n-1)){
		#cat(i)
		obslist[[i]] <- obj[[hord[i]]]$rlabs
		if (length(obslist[[i]])==2) {
			merge[i,] <- -obslist[[i]]
		} else {
			remaining <- obslist[[i]]
			pointer <- i
			idlist <- NULL
			while (length(remaining)>1){
				pointer <- pointer - 1
				#cat('\n',pointer,':',remaining,':',obslist[[pointer]],'\n')
				if (any(remaining %in% obslist[[pointer]])){
					idlist <- c(idlist,pointer)
					remaining <- remaining[!(remaining %in% obslist[[pointer]])]
				}
			}
			#browser()
			if (length(remaining)==1) idlist <- c(idlist,-remaining)
			if (length(idlist)!=2) {
				browser()
				stop('code error:wrong length list\n')
			}
			merge[i,] <- idlist
		}
		#cat('done i')
	}
	storage.mode(merge) <- 'integer'
  myorder <- descendants(merge,length(height))[[1]]
  storage.mode(myorder) <- 'integer'
	result <- list(height=height,merge=merge,order=myorder)
	class(result) <- 'hclust'
	result
}
