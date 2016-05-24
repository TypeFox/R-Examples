`boot.conf.values` <-
function(data, random.seed=12345, nrep=500) {
  dat <- as.matrix(data)
  naber <- ncol(dat)
  nrep <- nrep
  trees <- array(dim=c(2, 2*naber-1, nrep))
  ps <- matrix(NA, nrow=2*naber-1, ncol=nrep) #the parameters

  for (rep in 1:nrep) {
	set.seed(random.seed+rep)
#	cat("repetition", rep, "\n")
	permdat <- dat[sample(nrow(dat), replace=TRUE),]
	if (all(colSums(permdat) != 0)) {
		res <- MLtopology(permdat)
		trees[,,rep] <- res$tree
		ps[,rep] <- res$p
	}
  }

  splitlist <- matrix(nrow=0, ncol=naber+2)
  for (rep in 1:nrep) {
	tree <- trees[,,rep]
	for (edge in 1:ncol(tree)) { # record the split
		if (ps[edge,rep] < 1 & !is.na(ps[edge,rep])) {
			subt <- subtree(tree[2,edge], tree)
			if (length(subt)>0) {
				split <- rep(0, naber)
				split[which(1:naber %in% tree[2,subt])] <- 1
				int <- as.vector(split %*% 2^(0:(naber-1)))
				a <- as.vector(splitlist[,naber+1])
				splitlistind <- which(a==int)
				if (length(splitlistind)>0) {
					splitlist[splitlistind, naber+2] <- splitlist[splitlistind, naber+2]+1
				} else {
					splitlist <- rbind(splitlist, c(split, int, 1))
				}
			}
		}
	}
  }
  majsplits <- which(splitlist[,naber+2] > 0.1*nrep)
##display the splits occurring in > 10% of the bootstrap data sets (the splits are characterized by one of the two subsets of leaves):
  res.list <- list()
  res.list[["frequencies"]] <- splitlist[majsplits,]
  res.list[["confidence.values"]] <- splitlist[majsplits, ncol(splitlist)]/nrep
  names(res.list[["confidence.values"]]) <- apply(splitlist[majsplits,1:(ncol(splitlist)-2)], 1, function(x)paste(colnames(dat)[which(x==1)], collapse=" "))
  return(res.list)
}

