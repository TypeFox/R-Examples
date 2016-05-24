cat('
Demo for hybridHclust library
Commentary on the example is online at
 http://www.acadiau.ca/~hchipman/hybridHclust/tutorial.html')

readline('hit enter to start')

data(sorlie)
x <- scale(sorlie) 
data(sorlielabels)
dmat <- 1-cor(x)
dimnames(dmat) <- list(1:nrow(dmat),1:nrow(dmat))

# find and print MCs
mc1 <- mutualCluster(distances=dmat,plot=TRUE)
print(mc1)
get.distances(mc1)

# find hybrid hierarchical model and explore it
hyb1 <- hybridHclust(t(x),mc1,trace=TRUE)
plot(hyb1)
table(hyb=cutree(hyb1,5),sorlie=sorlielabels)

# get conventional bottom-up and top-down clusterings
h1 <- hclust(as.dist(dmat),method='average')
t1 <- tsvq(t(x),K=nrow(t(x)),row.labs=1:nrow(t(x)))

# comparative exploration
treelist <- vector('list',3)
treelist[[1]] <- list(tr=h1,name='Bottom-up')
treelist[[2]] <- list(tr=t1,name='Top-down')
treelist[[3]] <- list(tr=hyb1,name='Hybrid')

# define 4 helper functions....

make.tree.partitions <- function(treelist,prange){
# use: make.tree.partitions(list(list(tr=hclust(dist(x)),name='hclust')),2:20)
	parts <- vector('list',length(treelist))
	thenames <- NULL
	for (i in 1:length(treelist)){
  	parts[[i]] <- matrix(0,length(treelist[[i]]$tr$height)+1,max(prange))
  	for (j in prange) {
			if (j==1) parts[[i]][,j] <- rep(1,nrow(parts[[i]]))
			else if (j==nrow(parts[[i]])) parts[[i]][,j] <- 1:nrow(parts[[i]])
			else parts[[i]][,j] <- cutree(treelist[[i]]$tr,j)
		}
  	thenames <- c(thenames,treelist[[i]]$name)
	}
	list(parts,thenames,prange)
}
 
calc.partition.dists <- function(thelist,prange)
{
	parts <- thelist[[1]]
	thenames <- thelist[[2]]
	#prange <- thelist[[3]]
	np <- length(parts)
	dtreemat <- matrix(0,np,np)
	dimnames(dtreemat) <- list(thenames,thenames)
	tdists <- NULL
	for (ii in prange){
	for (i in 1:(np-1)) for (j in (i+1):np)
  	dtreemat[i,j] <- dtreemat[j,i] <-
    	dpart(parts[[i]][,ii],parts[[j]][,ii])
  	#print(as.dist(dtreemat))
  	tdists <- cbind(tdists,as.dist(dtreemat))
	}
	tnames <- NULL
	for (i in 1:(np-1)) for (j in (i+1):np)
  	tnames <- c(tnames,paste(thenames[i],thenames[j],sep=':'))
	rownames(tdists) <- tnames
	tdists
}

dpart <- function(p1,p2){
	# much more efficient version that uses cross-tab of node id's for two trees
	# p1 and p2 are vectors of cluster labels.
	tt <- table(p1,p2)
	d <- 0
	for (i in 1:nrow(tt))
		for (j in 1:ncol(tt)){
			if (i<nrow(tt)) d <- d+sum(tt[i,j]*tt[(i+1):nrow(tt),j])
			if (j<ncol(tt)) d <- d+sum(tt[i,j]*tt[i,(j+1):ncol(tt)])
		}
	return(d/choose(length(p1),2))
}

calc.sum.dist <- function(dmat,clus){
# note that clus is a matrix of indices, with each column a diff partition
  sum.d <- rep(0,ncol(clus))
  for (k in 1:ncol(clus)){
    actual <- rep(0,k)
    for (j in 1:k) {
      actual[j]<- sum(dmat[clus[,k]==j,clus[,k]==j])/2/sum(clus[,k]==j)
    }
    sum.d[k] <- sum(actual)
  }
  sum.d
}

# ...end of helper function definitions.

# Compare the resultant trees, via comparison of partitions, and
# within-group sum of distances
partitions <- make.tree.partitions(treelist,1:85)

whichplot <- 2:84
junk2 <- calc.partition.dists(partitions,whichplot)
par(mfrow=c(1,1))
matplot(whichplot,t(junk2),type='l',col=1,lwd=2,xlab='Number of clusters',
  ylab='dissimilarity between clusterings',ylim=c(0,max(junk2)),xlim=c(1,85)
  ,log='x')
matpoints(whichplot[1:9],t(junk2[,1:9]),pch=19,col=1)
legend(mean(whichplot)/10,max(junk2),dimnames(junk2)[[1]],
  lty=1:nrow(junk2),lwd=2)

tot.dist.from.mean <- matrix(0,length(partitions[[2]]),85)
for (i in 1:length(partitions[[2]])){
    tot.dist.from.mean[i,] <- calc.sum.dist(dmat,partitions[[1]][[i]])
}
par(mfrow=c(1,2))
mylty<-1:6
mycol <- c('black','red','blue')
matplot(t(tot.dist.from.mean),type='l',xlab='number of clusters',
  ylab='within cluster sum of distances',col=mycol,lty=mylty)
legend(20,max(tot.dist.from.mean),partitions[[2]],lty=mylty[1:3],col=mycol[1:3])
relative <- tot.dist.from.mean
for (i in 1:ncol(relative)) 
  relative[,i] <- relative[,i] - min(relative[,i])
matplot(t(relative),type='l',col=mycol,xlab='number of clusters',
  ylab='within cluster sum of distances (relative to min)',lty=mylty,lwd=1.5)
legend(20,max(relative),partitions[[2]],lty=mylty[1:3],col=mycol[1:3],lwd=1.5)

# comparison with sorlie labels
for (i in 1:3) 
  cat(treelist[[i]]$name, dpart(sorlielabels,cutree(treelist[[i]]$tr,5)),'\n')

