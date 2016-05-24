ConsensualWords <- function (res.fast, nbtimes = 3, nbsimul = 500, proba = 0.05, graph = TRUE, axes = c(1,2))
{
	if (nbtimes < 2) stop("nbtimes should be superior or equal to 2")
	ncp <- res.fast$call$ncp
	if (ncp < 2)     stop("ncp should be superior or equal to 2")

	mots <- rownames(res.fast$descriptor)
	coord.mot.acm <- res.fast$acm$var$coord
	nm <- tolower(rownames(coord.mot.acm))

	nb.times <- numeric()
	centroids <- matrix(0, ncol=ncp, 0)
	within.inertia <- numeric()

	WordinList <- function(mot,liste){ ## return a vector of boolean to know if the word belongs to the list 
	  aa=sapply(liste,strsplit,res.fast$call$sep.words)
#	  aa=sapply(liste,strsplit," ")
	  aa=lapply(lapply(aa,strsplit,"_"),unlist)
	  bb= sapply(lapply(aa,match,mot,nomatch=0),sum)
	  bb = as.logical(bb)
	}
	## Calculation of the barycenters and within inertias
	for (i in 1:length(mots)){
##		tab <- subset(coord.mot.acm, grepl(mots[i], nm))
		tab <- subset(coord.mot.acm, WordinList(mots[i], nm))
		nb.times[i] <- dim(tab)[1]
		bc <- apply(tab, 2, mean)
		centroids <- rbind(centroids,bc)
		within.inertia[i] <- round(sum(apply(tab, 2, var)),3)
	}

	rownames(centroids) <- mots
	within.inertia <- cbind(nb.times,within.inertia)
	rownames(within.inertia) <- mots

	## Checking the words according to their frequency of use
	filtrer <- rownames(within.inertia[nb.times >= nbtimes,])

	## Bootstrap
	prob <- numeric()
	for (i in 1:length(filtrer)){
		Rsamples <- replicate(nbsimul, sample(1:nrow(coord.mot.acm), within.inertia[filtrer[i],1], replace = TRUE))
		ii.bt <- numeric()
		for (j in 1:nbsimul){
			tab <- coord.mot.acm[Rsamples[,j],]
			ii.bt[j] <- sum(apply(tab, 2, var))
		}
		prob[i] <- mean(ii.bt <= within.inertia[filtrer[i],2])
	}

	mot.con <- cbind.data.frame(within.inertia[filtrer,], prob)
	rownames(mot.con) <- filtrer
	mot.con <- mot.con[order(mot.con[,3]),]
	Mots.consensuels <- rownames(mot.con[mot.con[,3] <= proba,])

	if (graph) {
     res.fast$acm$var$coord <- centroids[Mots.consensuels,]
	 plot(res.fast$acm,axes=axes,title="Sorting task with consensual words",new.plot=TRUE)
	}

	return(list(Centroids = centroids[filtrer,], Within.inertia = within.inertia[filtrer,], Results.Bootstrap = mot.con, Consensual.words = Mots.consensuels[order(Mots.consensuels)]))
}



