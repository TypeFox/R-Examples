seqlogp <- function(seqdata, prob="trate", time.varying=TRUE, begin="freq", weighted=TRUE) {
	## Liste des taux de transitions par age
	sl <- seqlength(seqdata)
	maxage <- max(sl)
	nbtrans <- maxage -1
	agedtr <- vector(mode="list", length=maxage)
	
	## On ajoute 1 pour que les codes correspondent aux index R (commence à 1)
	seqdatanum <- seqasnum(seqdata)+1
	nbstates <- max(seqdatanum)
	## User defined begin frequencies
	if(is.numeric(begin)){
		if (length(begin)!=nbstates) {
			stop("Begin frequencies should be a numeric vector of length ", nbstates)
		}
		message(" [>] Using user defined frequencies as starting point")
		firstfreq <- begin
	}
	##Compute from data
	else if (is.character(begin) && begin=="freq") {
		message(" [>] Using data defined frequencies as starting point")
		firstfreq <- seqstatd(seqdata, weighted=weighted)$Frequencies[, 1]
	}
	else {
		stop("Unknow method to compute starting frequencies")
	}
	
	###Automatic method to compute transition rates
	if (is.character(prob)) {
		if (prob=="trate") {
			if (time.varying) {
				message(" [>] Using time varying transition rates as probability model")
				agedtr <- seqtrate(seqdata, time.varying=TRUE, weighted=weighted)
			}
			else {
				message(" [>] Using global transition rates as probability model")
				agedtr <- array(0, dim=c(nbstates, nbstates, nbtrans))
				tr <- seqtrate(seqdata, weighted=weighted)
				for (i in 1:nbtrans) {
					agedtr[,,i] <- tr
				}
			}
		}
		else if (prob=="freq") {
			## On crée quand même une matrice de transition (qui ne dépend pas de l'état précédant)
			## On peut ainsi utiliser le même algorithme
			message(" [>] Using time varying frequencies as probability model")
			agedtr <- array(0, dim=c(nbstates, nbstates, nbtrans))
			if (time.varying) {
				freqs <- seqstatd(seqdata, weighted=weighted)$Frequencies
				for (i in 1:nbtrans) {
					for (j in 1:length(freqs[, i+1])) {
							agedtr[, j,i] <- freqs[j, i+1]
					}
				}
			}
			else {
				message(" [>] Using global frequencies as probability model")
				freqs <- seqstatf(seqdata, weighted=weighted)$Percent/100
				for (i in 1:nbtrans) {
					for (j in 1:length(freqs)) {
						agedtr[, j,i] <- freqs[j]
					}
				}
			}
		}
		else {
			stop("Unknow method to compute transition rate")
		}
	}
	## User defined transition rates
	else{
		if(is.array(prob)){
			if(length(dim(prob)) == 3) {
				##Correct dimensions
				if(any(dim(prob)!=c(nbstates, nbstates, nbtrans))){
					stop("Transition rate should be an array of size (state x state x transition) ",
						nbstates,"x", nbstates, "x", nbtrans)
				}
				message(" [>] Using user defined time varying transition rates as probability model")
				agedtr <- prob
			} else if (length(dim(prob)) == 2) {
				message(" [>] Using user defined global transition rates as probability model")
				if(any(dim(prob)!=c(nbstates, nbstates))){
					stop("Transition rate should be a matrix of size (state x state) ",
						nbstates,"x", nbstates)
				}
				agedtr <- array(0, dim=c(nbstates, nbstates, nbtrans))
				for (i in 1:nbtrans) {
					agedtr[,,i] <- prob
				}
			}
			else {
				stop("Transition rate should be an array of size (state x state x transition) ",
						nbstates,"x", nbstates, "x", nbtrans, " or a matrix of size (state x state) ",
						nbstates,"x", nbstates)
			}
		}
		else {
			stop("Unknow method to compute transition rate")
		}
	}
	logp <- numeric(length=(nrow(seqdata)))
	logp[] <- 0
	for (i in 1:nrow(seqdatanum)) {
		p <- firstfreq[seqdatanum[i, 1]]
		logp[[i]] <- -log(p)
		if (sl[i]>1) {
			for (j in 2:sl[i]) {
				p <- agedtr[seqdatanum[i, j-1], seqdatanum[i, j], j-1]
				logp[[i]] <- logp[[i]] -log(p)
			}
		}
	}
	return(logp)
}
