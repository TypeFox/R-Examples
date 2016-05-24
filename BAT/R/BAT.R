#####BAT - Biodiversity Assessment Tools
#####Version 1.5.0 (2016-04-21)
#####By Pedro Cardoso, Francois Rigal, Jose Carlos Carvalho
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P., Rigal, F. & Carvalho, J.C. (2015) BAT - Biodiversity Assessment Tools, an R package for the measurement and estimation of alpha and beta taxon, phylogenetic and functional diversity. Methods in Ecology and Evolution, 6, 232-236.
#####Changed from v1.4.0:
#####Added function optim.spatial
#####Added new SAR equations
#####Minor tweaks in some functions

#####BAT Stats:
#####library("cranlogs")
#####day <- cran_downloads(package = "BAT", from = "2014-08-19", to = "2015-10-19")
#####group <- matrix(day$count, 140, byrow=T)
#####plot(rowSums(group), type = "n")
#####lines(rowSums(group))

#####required packages
library("graphics")
library("nls2")
library("raster")
library("spatstat")
library("stats")
library("utils")
library("vegan")
#' @import graphics
#' @import nls2
#' @import spatstat
#' @import stats
#' @import utils
#' @import vegan
#' @importFrom raster rasterize
#' @importFrom raster rasterToPoints

#####auxiliary functions
prep <- function(comm, xtree, abund = TRUE){
	len <- xtree[[1]] 							## length of each branch
	A <- xtree[[2]]									## matrix species X branches
	minBranch <- min(len[colSums(A)==1]) 	## minimum branch length of terminal branches
	BA <- comm%*%A 												## matrix samples X branches
	if (!abund)	BA = ifelse(BA >= 1, 1, 0)
	return (list(lenBranch = len, sampleBranch = BA, speciesBranch = A, minBranch = minBranch))
}

rarefaction <- function(comm){
	n <- sum(comm)
	for (s in 1:nrow(comm))
		n <- min(n, sum(comm[s,]))
	return(n)
}

rss <- function(x, y){
	return (sum((x-y)^2))
}

AIC <- function(x, y, k){
	n = length(x)
	return(n * log(rss(x,y)/n) + 2*k)
}

AICc <- function(x, y, k){
	n = length(x)
	return(AIC(x, y, k) + (2*k*(k+1))/(n-k-1))
}

r2 <- function(x, y){
	SSn <- rss(x, y)
	SSd <- sum((y-mean(y))^2)
	return(1-(SSn/SSd))
}

logit <- function(x){
	return(log(x/(1-x)))
}

revLogit <- function(x){
	return(exp(x)/(1+exp(x)))
}

#####xTree function partly adapted from http://owenpetchey.staff.shef.ac.uk/Code/Code/calculatingfd_assets/Xtree.r
#####by Jens Schumacher (described in Petchey & Gaston 2002, 2006)
xTree <- function(tree) {
  if (class(tree) == "hclust"){
  	nSpp <- nrow(as.data.frame(tree['order']))
  	sppEdges <- matrix(0, nSpp, 2 * nSpp - 2) 
  	lenEdges <- vector("numeric", 2 * nSpp - 2)
  	for(i in 1:(nSpp - 1)) {
  		if(tree$merge[i, 1] < 0) {
  			lenEdges[2 * i - 1] <- tree$height[order(tree$height)[i]] 
  			sppEdges[ - tree$merge[i, 1], 2 * i - 1] <- 1
  		} else {
  			lenEdges[2 * i - 1] <- tree$height[order(tree$height)[i]] - tree$height[order(tree$height)[tree$merge[i, 1]]]
  			sppEdges[, 2 * i - 1] <- sppEdges[, 2 * tree$merge[i, 1] - 1] + sppEdges[ , 2 * tree$merge[i, 1]]
  		} 
  		if(tree$merge[i, 2] < 0) {
  			lenEdges[2 * i] <- tree$height[order(tree$height)[i]] 
  			sppEdges[ - tree$merge[i, 2], 2 * i] <- 1
  		} else {
  			lenEdges[2 * i] <- tree$height[order(tree$height)[i]] - tree$height[order(tree$height)[tree$merge[i, 2]]]
  			sppEdges[, 2 * i] <- sppEdges[, 2 * tree$merge[i, 2] - 1] + sppEdges[, 2 *tree$merge[i, 2]]
  		}
  	} 
  	rownames(sppEdges) <- tree$labels
  	list(lenEdges, sppEdges)
  } else if (class(tree) == "phylo"){
    lenEdges <- tree$edge.length
    nSpp <- length(tree$tip.label)
    nEdges <- length(tree$edge.length)
    root <- nSpp + 1
    sppEdges <- matrix(0, nSpp, nEdges)
    for(i in 1:nSpp){
      find = i                                    #start by finding the ith species
      repeat{
        row = which(tree$edge[,2] == find)        #locate in which row of the edge table is our species or edge to be found
        sppEdges[i, row] = 1
        find = tree$edge[row,1]                   #find next edge if any until reaching the root
        if(find == root) break                    #all edges of this species were found, go to next species
      }
    }
    rownames(sppEdges) <- tree$tip.label
    list(lenEdges, sppEdges)
  } else {
    cat("Unrecognized tree object!")
  }
}

#####observed diversity
sobs <- function(comm, xtree){
	#if (is.vector(comm))
		#comm = matrix(c(comm,rep(0,length(comm))),ncol=2)
	if (missing(xtree)){
		return(length(colSums(comm)[colSums(comm) > 0]))
	} else {
		data <- prep(comm, xtree)
		value <- ifelse (colSums(data$sampleBranch) > 0, 1, 0) # vector of observed branches
		return (sum(value*data$lenBranch))
	}
}

#####diversity of rare species for abundance - singletons, doubletons, tripletons, etc
srare <- function(comm, xtree, n = 1){
	if(missing(xtree)){
		return(length(colSums(comm)[colSums(comm) == n]))
	} else {
		data <- prep(comm, xtree)
		value <- ifelse (colSums(data$sampleBranch) == n, 1, 0) # vector of branches with given abundance
		return (sum(value*data$lenBranch))
	}
}

#####diversity of rare species for incidence - uniques, duplicates, triplicates, etc
qrare <- function(comm, xtree, n = 1){
	if(missing(xtree)){
		comm <- ifelse(comm > 0, 1, 0)
		return(length(colSums(comm)[colSums(comm) == n]))
	} else {
		data <- prep(comm, xtree, FALSE)
		value <- ifelse (colSums(data$sampleBranch) == n, 1, 0) # vector of branches with given incidence
		return (sum(value*data$lenBranch))
	}
}

#####minimum terminal branch length, = 1 in case of TD
minBranch <- function(comm, xtree){
	if (missing(xtree)){
		return(1)
	} else {
		data <- prep(comm, xtree)
		return(data$minBranch)
	}
}

#####non-parametric estimators
chao <- function(obs, s1, s2, mb){
	return(obs + (s1*(s1-mb))/(2*(s2+mb)))
}

jack1ab <- function(obs, s1){
	return(obs + s1)
}

jack1in <- function(obs, q1, q){
	return(obs + q1 * ((q-1)/q))
}

jack2ab <- function(obs, s1, s2){
	return(obs + 2*s1 - s2)
}

jack2in <- function(obs, q1, q2, q){
	if (q > 1)	return(obs + (q1*(2*q-3)/q - q2*(q-2)^2/(q*(q-1))))
	else return(obs + 2*q1 - q2)
}

pcorr <- function(obs, s1){
	return(1+(s1/obs)^2)
}

#####observed beta (a = shared species/edges, b/c = species/edges exclusive to either site, comm is a 2sites x species matrix)
betaObs <- function(comm, xtree, abund = FALSE, func = "jaccard"){
  if(sum(comm) == 0)                                ##if no species on any community return 0
    return(list(Btotal = 0, Brepl = 0, Brich = 0))
	if (!abund || max(comm) == 1) {										##if incidence data
		obs1 <- sobs(comm[1,,drop=FALSE], xtree)
		obs2 <- sobs(comm[2,,drop=FALSE], xtree)
		obsBoth <- sobs(comm, xtree)
		a <- obs1 + obs2 - obsBoth
		b <- obsBoth - obs2
		c <- obsBoth - obs1
	} else if (abund & missing(xtree)){								##if abundance data
		a <- 0
		b <- 0
		c <- 0
		for (i in 1:ncol(comm)){
		  minComm <- min(comm[1,i], comm[2,i])
		  a <- a + minComm
		  b <- b + comm[1,i] - minComm
		  c <- c + comm[2,i] - minComm
		}
	} else {																					##if abundance and tree
		##due to the way Soerensen doubles the weight of the a component, using a tree or not will be the same with abundance data.
		data <- prep(comm, xtree)
		a = sum(data$lenBranch * apply(data$sampleBranch,2,min))
		diff = data$lenBranch * (data$sampleBranch[1,] - data$sampleBranch[2,])
		b = sum(replace(diff, diff < 0, 0))
		c = sum(replace(diff, diff > 0, 0) * -1)
	}
	denominator <- a + b + c
	if(tolower(substr(func, 1, 1)) == "s")
		denominator <- denominator + a
	return(list(Btotal = (b+c)/denominator, Brepl = 2*min(b,c)/denominator, Brich = abs(b-c)/denominator))
}

##create latitude layer
raster.lat <- function(layers){
	lat <- layers[[1]]
	x <- rasterToPoints(lat)[,1:2]
	lat <- rasterize(x[,1:2], lat, x[,2])
	names(lat) <- "latitude"
	return(lat)
}

##create longitude layer
raster.long <- function(layers){
	long <- layers[[1]]
	x <- rasterToPoints(long)[,1:2]
	long <- rasterize(x[,1:2], long, x[,1])
	names(long) <- "longitude"
	return(long)
}

##################################################################################
##################################MAIN FUNCTIONS##################################
##################################################################################


#' Alpha diversity (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).
#' @description Observed alpha diversity with possible rarefaction, multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details TD is equivalent to species richness. Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' The rarefaction option is useful to compare communities with much different numbers of individuals sampled, which might bias diversity comparisons (Gotelli & Colwell 2001)
#' @return A matrix of sites x diversity values (either "Obs" OR "Avg, Min, LowerCL, UpperCL and Max").
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Gotelli, N.J. & Colwell, R.K. (2001) Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecology Letters, 4, 379-391.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @examples comm <- matrix(c(0,0,1,1,0,0,2,1,0,0), nrow = 2, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha(comm)
#' alpha(comm, raref = 0)
#' alpha(comm, tree, 2, 100)
#' @export
alpha <- function(comm, tree, raref = 0, runs = 100){
  
  comm <- as.matrix(comm)
  if (!missing(tree))
		tree <- xTree(tree)
	
  nComm <- nrow(comm)
	if(raref < 1){						# no rarefaction if 0 or negative
		results <- matrix(0, nComm, 1)
		for (s in 1:nComm){
			results[s,1] <- sobs(comm[s,, drop=FALSE], tree)
		}
		rownames(results) <- rownames(comm)
		colnames(results) <- "Obs"
		return (results)
	}
	if (raref == 1)
		raref <- rarefaction(comm)				# rarefy by minimum n among all communities
	results <- matrix(0, nComm, 5)
	for (s in 1:nComm){
		res <- c()
		for (r in 1:runs){
			res <- c(res,sobs(rrarefy(comm[s,], raref), tree))
		}
		results[s,] <- c(mean(res), min(res), quantile(res, 0.025), quantile(res, 0.975), max(res))
	}
	rownames(results) <- rownames(comm)
	colnames(results) <- c("Avg", "Min", "LowerCL", "UpperCL", "Max")
	return (results)
}

#' Alpha diversity accumulation curves (observed and estimated).
#' @description Estimation of alpha diversity of a single site with accumulation of sampling units.
#' @param comm A sampling units x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for Phylogenetic (PD) or Functional (FD) Diversity, not for Taxon Diversity (TD)).
#' @param func The class of estimators to be used:
#' If func is partial match of "curve", TD, PD or FD are based on extrapolating the accumulation curve of observed diversity.
#' If func is partial match of "nonparametric", TD, PD or FD are based on non-parametric estimators.
#' If func is partial match of "completeness", PD or FD estimates are based on the completeness of TD (requires a tree to be used).
#' If not specified, default is "nonparametric.
#' @param target True diversity value to calculate the accuracy of curves (scaled mean squared error). If not specified do not calculate accuracy (default), -1 uses the total observed diversity as true diversity and any other value is the true known diversity.
#' @param runs Number of random permutations to be made to the sampling order. If not specified, default is 100.
#' @param prog Present a text progress bar in the R console.
#' @details Observed diversity often is an underestimation of true diversity. Several approaches have been devised to estimate species richness (TD) from incomplete sampling.
#' These include: (1) fitting asymptotic functions to randomised accumulation curves (Soberon & Llorente 1993; Flather 1996; Cardoso et al. in prep.)
#' (2) the use of non-parametric estimators based on the incidence or abundance of rare species (Heltshe & Forrester 1983; Chao 1984, 1987; Colwell & Coddington 1994).
#' A correction to non-parametric estimators has also been recently proposed, based on the proportion of singleton or unique species
#' (species represented by a single individual or in a single sampling unit respectively; Lopez et al. 2012).
#' Cardoso et al. (2014) have proposed a way of adapting these approaches to estimate PD and FD, also adding a third possible approach for
#' these dimensions of diversity: (3) correct PD and FD values based on the completeness of TD, where completeness equals the proportion of estimated true diversity that was observed.
#' Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of sampling units x diversity values (sampling units, individuals, observed and estimated diversity).
#' The values provided by this function are: 
#' @return Sampl - Number of sampling units;
#' @return Ind - Number of individuals;
#' @return Obs - Observed diversity;
#' @return S1 - Singletons;
#' @return S2 - Doubletons;
#' @return Q1 - Uniques;
#' @return Q2 - Duplicates;
#' @return Jack1ab - First order jackknife estimator for abundance data;
#' @return Jack1in - First order jackknife estimator for incidence data;
#' @return Jack2ab - Second order jackknife estimator for abundance data;
#' @return Jack2in - Second order jackknife estimator for incidence data;
#' @return Chao1 - Chao estimator for abundance data;
#' @return Chao2 - Chao estimator for incidence data;
#' @return Clench - Clench or Michaelis-Menten curve;
#' @return Exponential - Exponential curve;
#' @return Rational - Rational function;
#' @return Weibull - Weibull curve;
#' @return The P-corrected version of all non-parametric estimators is also provided.
#' @return Accuracy - if accuracy is to be calculated a list is returned instead, with the second element being the scaled mean squared error of each estimator.
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, in press.
#' @references Chao, A. (1984) Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.
#' @references Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' @references Colwell, R.K. & Coddington, J.A. (1994). Estimating terrestrial biodiversity through extrapolation. Phil. Trans. Roy. Soc. London B 345, 101-118.
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Flather, C. (1996) Fitting species-accumulation functions and assessing regional land use impacts on avian diversity. Journal of Biogeography, 23, 155-168.
#' @references Heltshe, J. & Forrester, N.E. (1983) Estimating species richness using the jackknife procedure. Biometrics, 39, 1-11.
#' @references Lopez, L.C.S., Fracasso, M.P.A., Mesquita, D.O., Palma, A.R.T. & Riul, P. (2012) The relationship between percentage of singletons and sampling effort: a new approach to reduce the bias of richness estimates. Ecological Indicators, 14, 164-169.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @references Soberon, M.J. & Llorente, J. (1993) The use of species accumulation functions for the prediction of species richness. Conservation Biology, 7, 480-488.
#' @examples comm <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha.accum(comm)
#' alpha.accum(comm, func = "nonparametric")
#' alpha.accum(comm, tree, "completeness")
#' alpha.accum(comm, tree, "curve", runs = 1000)
#' alpha.accum(comm, target = -1)
#' @export
alpha.accum <- function(comm, tree, func = "nonparametric", target = -2, runs = 100, prog = TRUE){
	
	comm <- as.matrix(comm)
	if (!missing(tree))
		tree <- xTree(tree)
	
	#####function options:
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	#####completeness (PD/FD with TD completeness correction)
	#####curve (TD/PD/FD with curve fitting)
	func <- match.arg(func, c("nonparametric", "completeness", "curve"))
	
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	switch(func, nonparametric = {
		resultsArray <- array(0, dim = c(nrow(comm), 19, runs))
		if(target > -2)
			smse <- matrix(0, runs, 19)
		if (prog) pb <- txtProgressBar(0, runs, style = 3)
		for (r in 1:runs){
			comm <- comm[sample(nrow(comm)),, drop=FALSE]			#shuffle rows (sampling units)
			data <- matrix(0,1,ncol(comm))
			runData <- matrix(0,nrow(comm),19)
			for (q in 1:nrow(comm)){
				data <- rbind(data, comm[q,])
				n <- sum(rowSums(data))
				obs <- sobs(data, tree)
				s1 <- srare(data, tree, 1)
				s2 <- srare(data, tree, 2)
				q1 <- qrare(data, tree, 1)
				q2 <- qrare(data, tree, 2)
				mb <- minBranch(data, tree)
				j1ab <- jack1ab(obs, s1)
				j1abP <- j1ab * pcorr(obs, s1)
				j1in <- jack1in(obs, q1, q)
				j1inP <- j1in * pcorr(obs, q1)
				j2ab <- jack2ab(obs, s1, s2)
				j2abP <- j2ab * pcorr(obs, s1)
				j2in <- jack2in(obs, q1, q2, q)
				j2inP <- j2in * pcorr(obs, q1)
				c1 <- chao(obs, s1, s2, mb)
				c1P <- c1 * pcorr(obs, s1)
				c2 <- chao(obs, q1, q2, mb)
				c2P <- c2 * pcorr(obs, q1)
				runData[q,] <- c(q, n, obs, s1, s2, q1, q2, j1ab, j1abP, j1in, j1inP, j2ab, j2abP, j2in, j2inP, c1, c1P, c2, c2P)
			}
			resultsArray[,,r] <- runData
			if(exists("smse")){					##if accuracy is to be calculated
				if(r == 1){
					if(target == -1){
						truediv <- runData[nrow(runData),3]
					}else{
						truediv <- target
					}
				}
				s <- accuracy(runData, truediv)
				smse[r,3] <- s[1,1]
				smse[r,8:19] <- s[1,-1]
			}
			if (prog) setTxtProgressBar(pb, r)
		}
		if (prog) close(pb)
		
		#####calculate averages or medians of all runs
		results <- matrix(0,nrow(comm),19)
		v <- array(0, dim = c(runs))
		for (i in 1:nrow(comm)){
			for (j in 1:19){
				for (k in 1:runs){
					v[k] <- resultsArray[i,j,k]
				}
				if (j < 16 || missing(tree))
					results[i,j] <- mean(v)
				else
					results[i,j] <- median(v)
			}
		}
		if(exists("smse"))						##calculate accuracy
			smse <- colMeans(smse)
		
		#####completeness (PD/FD with TD completeness correction)
	}, completeness = {
		if (missing(tree))
			stop("Completeness option not available without a tree...")
		results <- alpha.accum(comm, runs = runs)
		obs <- matrix(0,nrow(comm),1)
		for (r in 1:runs){
			comm <- comm[sample(nrow(comm)),, drop=FALSE]			#shuffle rows (sampling units)
			for (s in 1:nrow(comm)){
				obs[s,1] <- obs[s,1] + sobs(comm[1:s,], tree)
			}
		}
		obs <- obs / runs
		for (i in 8:19)
			results[,i] <- obs * (results[,i] / results[,3])
		results[,3] <- obs
		
		#####curve (TD/PD/FD with curve fitting)
	}, curve = {
    results <- matrix(NA,nrow(comm),7)
    results[,1] <- seq(1,nrow(comm))  ##fill samples column
    results[,2] <- seq(sum(comm)/nrow(comm),sum(comm), sum(comm)/nrow(comm))  ##fill individuals column
    runObs <- rep(0,nrow(comm))
    if (prog) pb <- txtProgressBar(0, runs, style = 3)
    for (r in 1:runs){
 		  comm <- comm[sample(nrow(comm)),, drop=FALSE]		#shuffle rows (sampling units)
 		  for (s in 1:nrow(comm)){
 		    runObs[s] <- runObs[s] + sobs(comm[1:s,,drop=FALSE], tree)
 		  }
 		  if (prog) setTxtProgressBar(pb, r)
    }
    if (prog) close(pb)
		results[,3] <- runObs / runs
		
    rich <- results[nrow(comm),3]

		for (s in 3:nrow(results)){				##fit curves only with 3 or more sampling units
		  ## curve fitting
      x <- results[1:s,1]
			y <- results[1:s,3]
			##Clench
			stlist <- data.frame(a = rich, b = c(0.1, 0.5, 1))
			form <- y ~ (a*x)/(b+x)
			mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
			if(class(curve) != "try-error"){
				a <- coef(curve)[1]
				results[s,5] <- a
			}
			##Negative exponential
      form <- y ~ a*(1-exp(-b*x))
			mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
			if(class(curve) != "try-error"){
				a <- coef(curve)[1]
				results[s,6] <- a
			}
			##Rational
			stlist <- data.frame(a = rich, b = c(0.1, 0.5, 1, 5, 10), c = c(1, 10, 100, 1000, 10000))
			form <- y ~ (c+(a*x))/(b+x)
			mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
			if(class(curve) != "try-error"){
				a <- coef(curve)[1]
				results[s,4] <- a
			}
			##Weibull
 			stlist <- data.frame(a = rich, b = c(0,1,10), c = c(0,0.1,1))
 			form <- y ~ a*(1-exp(-b*(x^c)))
      mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
 			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
  			if(class(curve) != "try-error"){
  				a <- coef(curve)[1]
  				results[s,7] <- a
  			}
		}
		colnames(results) <- c("Sampl", "Ind", "Obs", "Clench", "Exponential", "Rational", "Weibull")
		return (results)
	})
	colnames(results) <- c("Sampl", "Ind", "Obs", "S1", "S2", "Q1", "Q2", "Jack1ab", "Jack1abP", "Jack1in", "Jack1inP", "Jack2ab", "Jack2abP", "Jack2in", "Jack2inP", "Chao1", "Chao1P", "Chao2", "Chao2P")
	if(exists("smse")){
		smse <- matrix(smse,1,length(smse))
		colnames(smse) <- colnames(results)
		smse[,c(1:2,4:7)] <- NA
		return(list(results, smse))
	}	else {
		return(results)
	}
}

#' Alpha diversity estimates.
#' @description Estimation of alpha diversity of multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundances or number of incidences.
#' @param tree An hclust or phylo object (used only for Phylogenetic (PD) or Functional (FD) Diversity, not for Taxon Diversity (TD)).
#' @param func The class of estimators to be used:
#' If func is partial match of "nonparametric", TD, PD or FD are based on non-parametric estimators.
#' If func is partial match of "completeness", PD or FD estimates are based on the completeness of TD (requires a tree to be used).
#' If not specified, default is "nonparametric".
#' @details Observed diversity often is an underestimation of true diversity.
#' Non-parametric estimators based on the incidence or abundance of rare species have been proposed to overcome the problem of undersampling (Heltshe & Forrester 1983; Chao 1984, 1987; Colwell & Coddington 1994).
#' A correction to non-parametric estimators has also been recently proposed, based on the proportion (P) of singleton or unique species
#' (species represented by a single individual or in a single sampling unit respectively; Lopez et al. 2012).
#' Cardoso et al. (2014) have proposed a way of adapting non-parametric species richness estimators to PD and FD. They have also proposed correcting PD and FD values based on the completeness of TD, where completeness equals the proportion of estimated true diversity that was observed.
#' Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of sites x diversity values (individuals, observed and estimated diversity).
#' The values provided by this function are: 
#' @return Ind - Number of individuals;
#' @return Obs - Observed diversity;
#' @return S1 - Singletons;
#' @return S2 - Doubletons;
#' @return Jack1ab - First order jackknife estimator for abundance data;
#' @return Jack2ab - Second order jackknife estimator for abundance data;
#' @return Chao1 - Chao estimator for abundance data.
#' @return The P-corrected version of all estimators is also provided.
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, in press.
#' @references Chao, A. (1984) Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.
#' @references Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' @references Colwell, R.K. & Coddington, J.A. (1994). Estimating terrestrial biodiversity through extrapolation. Phil. Trans. Roy. Soc. London B 345, 101-118.
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Heltshe, J. & Forrester, N.E. (1983) Estimating species richness using the jackknife procedure. Biometrics, 39, 1-11.
#' @references Lopez, L.C.S., Fracasso, M.P.A., Mesquita, D.O., Palma, A.R.T. & Riul, P. (2012) The relationship between percentage of singletons and sampling effort: a new approach to reduce the bias of richness estimates. Ecological Indicators, 14, 164-169.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @examples comm <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha.estimate(comm)
#' alpha.estimate(comm, tree)
#' alpha.estimate(comm, tree, func = "completeness")
#' @export
alpha.estimate <- function(comm, tree, func = "nonparametric"){
	
  comm <- as.matrix(comm)
  if (max(comm) == 1)
		stop("No estimates are possible without abundance or incidence frequency data")
	if (!missing(tree))
		tree <- xTree(tree)
	
	#####function options:
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	#####completeness (PD/FD with TD completeness correction)
	func <- match.arg(func, c("nonparametric", "completeness"))
	
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	switch(func, nonparametric = {
		results <- matrix(0,0,10)
		for (s in 1:nrow(comm)){
			data <- comm[s,,drop = FALSE]
			obs <- sobs(data, tree)
			n <- sum(data)
			s1 <- srare(data, tree, 1)
			s2 <- srare(data, tree, 2)
			mb <- minBranch(data, tree)
			j1ab <- jack1ab(obs, s1)
			j1abP <- j1ab * pcorr(obs, s1)
			j2ab <- jack2ab(obs, s1, s2)
			j2abP <- j2ab * pcorr(obs, s1)
			c1 <- chao(obs, s1, s2, mb)
			c1P <- c1 * pcorr(obs, s1)
			results <- rbind(results, c(n, obs, s1, s2, j1ab, j1abP, j2ab, j2abP, c1, c1P))
		}
		
		#####completeness (PD/FD with TD completeness correction)
	}, completeness = {
		if (missing(tree))
			stop("Completeness option not available without a tree...")
		results <- alpha.estimate(comm, , "nonparametric")
		obs <- matrix(0,nrow(comm),1)
		for (s in 1:nrow(comm))
			obs[s,1] <- obs[s,1] + sobs(comm[s,], tree)
		for (i in 5:10)
			results[,i] <- obs[,1] * (results[,i] / results[,2])
		results[,2] <- obs[,1]
	})
	rownames(results) <- rownames(comm)
	colnames(results) <- c("Ind", "Obs", "S1", "S2", "Jack1ab", "Jack1abP", "Jack2ab", "Jack2abP", "Chao1", "Chao1P")
	return(results)
}

#' Beta diversity (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).
#' @description Beta diversity with possible rarefaction, multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for PD or FD).
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis. If not specified, default is FALSE.
#' @param func Partial match indicating whether the Jaccard or Soerensen family of beta diversity measures should be used.  If not specified, default is Jaccard.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich.
#' Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone; Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' The rarefaction option is useful to compare communities with much different numbers of individuals sampled, which might bias diversity comparisons (Gotelli & Colwell 2001)
#' @return Three distance matrices between sites, one per each of the three beta diversity measures (either "Obs" OR "Avg, Min, LowerCL, UpperCL and Max").
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Gotelli, N.J. & Colwell, R.K. (2001) Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecology Letters, 4, 379-391.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
#' @examples comm <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta(comm)
#' beta(comm, func = "Soerensen")
#' beta(comm, tree)
#' beta(comm, raref = 1)
#' beta(comm, tree, abund = TRUE, "s", raref = 2)
#' @export
beta <- function(comm, tree, abund = FALSE, func = "jaccard", raref = 0, runs = 100){
	
  comm <- as.matrix(comm)
  if (!missing(tree))
		tree <- xTree(tree)
	nComm <- nrow(comm)
	
	if(raref < 1){						# no rarefaction if 0 or negative
		results <- array(0, dim=c(nComm, nComm, 3))
		for (i in 1:(nComm-1)){
			for (j in (i+1):nComm){
				commBoth <- as.matrix(rbind(comm[i,], comm[j,]))
				betaValues <- betaObs(commBoth, tree, abund, func)
				results[j,i,] <- unlist(betaValues)
			}
		}
		results <- list(Btotal = as.dist(results[,,1]),Brepl = as.dist(results[,,2]),Brich = as.dist(results[,,3]))
		return (results)
	}
	if (raref == 1)
		raref <- rarefaction(comm)				# rarefy by minimum n among all communities
	results <- array(0, dim=c(nComm, nComm, 3, 5))
	
	for (i in 1:(nComm-1)){
		for (j in (i+1):nComm){
			run <- matrix(0, runs, 3)
			for (r in 1:runs){
				commBoth <- as.matrix(rbind(rrarefy(comm[i,], raref), rrarefy(comm[j,], raref)))
				betaValues <- betaObs(commBoth, tree, abund, func)
				run[r,1] <- betaValues$Btotal
				run[r,2] <- betaValues$Brepl
				run[r,3] <- betaValues$Brich
			}
			for (b in 1:3){
				results[j,i,b,1] <- mean(run[,b])
				results[j,i,b,2] <- min(run[,b])
				results[j,i,b,3] <- quantile(run[,b], 0.025)
				results[j,i,b,4] <- quantile(run[,b], 0.975)
				results[j,i,b,5] <- max(run[,b])
			}
		}
	}
	results.total <- list(Btotal = as.dist(results[,,1,1]), Btotal.min = as.dist(results[,,1,2]), Btotal.lowCL = as.dist(results[,,1,3]), Btotal.upCL = as.dist(results[,,1,4]), Btotal.max = as.dist(results[,,1,5]))
	results.repl <- list(Brepl = as.dist(results[,,2,1]), Brepl.min = as.dist(results[,,2,2]), Brepl.lowCL = as.dist(results[,,2,3]), Brepl.upCL = as.dist(results[,,2,4]), Brepl.max = as.dist(results[,,2,5]))
	results.rich <- list(Brich = as.dist(results[,,3,1]), Brich.min = as.dist(results[,,3,2]), Brich.lowCL = as.dist(results[,,3,3]), Brich.upCL = as.dist(results[,,3,4]), Brich.max = as.dist(results[,,3,5]))
	results <- c(results.total, results.repl, results.rich)
	return (results)
}

#' Beta diversity accumulation curves.
#' @description Beta diversity between two sites with accumulation of sampling units.
#' @param comm1 A sampling units x species matrix for the first site, with either abundance or incidence data.
#' @param comm2 A sampling units x species matrix for the second site, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for Phylogenetic (PD) or Functional (FD) Diversity, not for Taxon Diversity (TD)).
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis. If not specified, default is FALSE.
#' @param func Partial match indicating whether the Jaccard or Soerensen family of beta diversity measures should be used. If not specified, default is jaccard.
#' @param runs Number of random permutations to be made to the sampling order. If not specified, default is 100.
#' @param prog Present a text progress bar in the R console.
#' @details As widely recognized for species richness, beta diversity is also biased when communities are undersampled.
#' Beta diversity accumulation curves have been proposed by Cardoso et al. (2009) to test if beta diversity has approached an asymptote when comparing two undersampled sites.
#' The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich.
#' Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone;
#' Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm1 and comm2 must be the same as in tree. Also, the number of sampling units should be similar in both sites.
#' @return Three matrices of sampling units x diversity values, one per each of the three beta diversity measures (sampling units, individuals and observed diversity).
#' @references Cardoso, P., Borges, P.A.V. & Veech, J.A. (2009) Testing the performance of beta diversity measures based on incidence data: the robustness to undersampling. Diversity and Distributions, 15, 1081-1090.
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta.accum(comm1, comm2)
#' beta.accum(comm1, comm2, func = "Soerensen")
#' beta.accum(comm1, comm2, tree)
#' beta.accum(comm1, comm2, abund = TRUE)
#' beta.accum(comm1, comm2, tree, TRUE)
#' @export
beta.accum <- function(comm1, comm2, tree, abund = FALSE, func = "jaccard", runs = 100, prog = TRUE){
	
  if(nrow(comm1) < 2 || nrow(comm1) != nrow(comm2))
		stop("Both communities should have multiple and the same number of sampling units")
  comm1 <- as.matrix(comm1)
  comm2 <- as.matrix(comm2)
  if (!missing(tree))
		tree <- xTree(tree)
	
	nSamples <- nrow(comm1)
	results <- matrix(0,nSamples, 4)
	colnames(results) <- c("Sampl", "Btotal", "Brepl", "Brich")
	if (prog) pb <- txtProgressBar(0, runs, style = 3)
	for (r in 1:runs){
		comm1 <- comm1[sample(nSamples),, drop=FALSE]			#shuffle sampling units of first community
		comm2 <- comm2[sample(nSamples),, drop=FALSE]			#shuffle sampling units of second community
		for (q in 1:nSamples){
			commBoth <- as.matrix(rbind(colSums(comm1[1:q,,drop=FALSE]),colSums(comm2[1:q,,drop=FALSE])))
			results[q,1] <- results[q,1] + q
			betaValues <- betaObs(commBoth, tree, abund, func)
			results[q,2] <- results[q,2] + betaValues$Btotal
			results[q,3] <- results[q,3] + betaValues$Brepl
			results[q,4] <- results[q,4] + betaValues$Brich
		}
		if (prog) setTxtProgressBar(pb, r)
	}
	if (prog) close(pb)
	results <- results/runs
	return(results)
}

#' Beta diversity among multiple communities.
#' @description Beta diversity with possible rarefaction - multiple sites measure calculated as the average or variance of all pairwise values.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only for Phylogenetic (PD) or Functional (FD) Diversity, not for Taxon Diversity (TD)).
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.  If not specified, default is FALSE.
#' @param func Indicates whether the Jaccard or Soerensen family of beta diversity measures should be used. If not specified, default is jaccard.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details Beta diversity of multiple sites simultaneously is calculated as either the average or the variance among all pairwise comparisons (Legendre, 2014).
#' The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich.
#' Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone;
#' Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of beta measures x diversity values (average and variance).
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Legendre, P. (2014) Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography, in press.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
#' @examples comm <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta.multi(comm)
#' beta.multi(comm, func = "Soerensen")
#' beta.multi(comm, tree)
#' beta.multi(comm, raref = 1)
#' beta.multi(comm, tree, TRUE, "s", raref = 2)
#' @export
beta.multi <- function(comm, tree, abund = FALSE, func = "jaccard", raref = 0, runs = 100){
	pairwise <- beta(comm, tree, abund, func, raref, runs)
	Btotal.avg <- mean(pairwise$Btotal)
	Brepl.avg <- mean(pairwise$Brepl)
	Brich.avg <- mean(pairwise$Brich)
	Btotal.var <- sum(pairwise$Btotal)/(ncol(comm)*(ncol(comm)-1))
	Brepl.var <- sum(pairwise$Brepl)/(ncol(comm)*(ncol(comm)-1))
	Brich.var <- sum(pairwise$Brich)/(ncol(comm)*(ncol(comm)-1))
	results <- matrix(c(Btotal.avg, Brepl.avg, Brich.avg, Btotal.var, Brepl.var, Brich.var), nrow = 3, ncol = 2)
	colnames(results) <- c("Average", "Variance")
	rownames(results) <- c("Btotal", "Brepl", "Brich")
	return(results)
}

#' Contribution of individuals or species to total PD or FD.
#' @description Contribution of each individual or species to the total PD or FD of a number of communities.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object.
#' @param abund A boolean (T/F) indicating whether contribution should be calculated per individual (T) or species (F). If not specified, default is TRUE.
#' @param relative A boolean (T/F) indicating whether contribution should be relative to total PD or FD (proportional contribution per individual or species). If False, the sum of contributions for each site is equal to total PD/FD, if True it is 1.
#' @details Contribution is equivalent to the evolutionary distinctiveness index (ED) of Isaac et al. (2007) if done by species and to the abundance weighted evolutionary distinctiveness (AED) of Cadotte et al. (2010) if done by individual.
#' @return A matrix of sites x species values.
#' @references Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. & Baillie, J.E.M. (2007) Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS One, 2: e296.
#' @references Cadotte, M.W., Davies, T.J., Regetz, J., Kembel, S.W., Cleland, E. & Oakley, T.H. (2010) Phylogenetic diversity metrics for ecological communities: integrating species richness, abundance and evolutionary history. Ecology Letters, 13: 96-105.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' contribution(comm, tree)
#' contribution(comm, tree, FALSE)
#' contribution(comm, tree, abund = FALSE, relative = TRUE)
#' @export
contribution <- function(comm, tree, abund = TRUE, relative = FALSE){
  if(!abund)
		comm <- ifelse(comm > 0, 1, 0)
  if(missing(tree))
    tree = hclust(as.dist(matrix(1,ncol(comm),ncol(comm))))
  if(class(tree) == "hclust")
    nEdges <- length(tree$merge)
  else
    nEdges <- length(tree$edge.length)
  comm <- as.matrix(comm)
	contrib <- matrix(0,nrow(comm),ncol(comm))
  
  for (i in 1:nrow(comm)){											#cycle through all sites/samples
		dataSample <- prep(comm[i,], xTree(tree), TRUE)
		valueBranch <- dataSample$lenBranch / dataSample$sampleBranch
		valueBranch <- ifelse(valueBranch == Inf, 0, valueBranch)
		for (j in 1:ncol(comm)){										#cycle through all species
			for (k in 1:nEdges){	            				#cycle through all branches
				contrib[i,j] <- contrib[i,j] + dataSample$speciesBranch[j,k] * valueBranch[k] * comm[i,j]
			}
		}
	}
  if(abund){
    contrib <- contrib/comm                      #so that the contribution per individual is given if it is the case
    contrib <- ifelse(contrib == "NaN", 0, contrib)
  }
  if(relative)
    contrib <- apply(contrib, 2, function(x) x / alpha(comm,tree))
	return(contrib)
}

#' Phylogenetic/functional dispersion of individuals or species.
#' @description Average dissimilarity between any two individuals or species randomly chosen in a community with replacement.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object.
#' @param abund A boolean (T/F) indicating whether dissimilarity should be calculated per individual (T) or species (F). If not specified, default is TRUE.
#' @param relative A boolean (T/F) indicating whether dissimilarity should be relative to the maximum distance between any two species in the tree.
#' @details If abundance data is used and a tree is given, dispersion is the quadratic entropy of Rao (1982).
#' If abundance data is not used but a tree is given, dispersion is the phylogenetic dispersion measure of Webb et al. (2002) although with replacement.
#' If abundance data is used but no tree is given, dispersion is 1 - Simpson's index (Simpson 1949).
#' @return A vector of values per site.
#' @references Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21: 24-43.
#' @references Simpson, E.H. (1949) Measurement of diversity. Nature 163: 688.
#' @references Webb, C.O., Ackerly, D.D., McPeek, M.A. & Donoghue, M.J. (2002) Phylogenies and community ecology. Annual Review of Ecology and Systematics, 33: 475-505.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' dispersion(comm)
#' dispersion(comm, tree)
#' dispersion(comm, tree, abund = FALSE)
#' dispersion(comm, tree, abund = FALSE, relative = TRUE)
#' @export
dispersion <- function(comm, tree, abund = TRUE, relative = FALSE){
	if(!abund)
		comm <- ifelse(comm > 0, 1, 0)
	if(missing(tree))
		tree = hclust(as.dist(matrix(1,ncol(comm),ncol(comm))))
	disp <- rep(0,nrow(comm))

  unique <- uniqueness(comm, tree, abund, relative)
  for (i in 1:nrow(comm)){  						                         #cycle through all sites/samples
  	present <- which(comm[i,]>0)                                 #which species exist in this site
  	proportion <- comm[i,present]/sum(comm[i,])                  #proportion incidence/abundance of species in this site
    disp[i] <- sum(unique[i,present]*proportion)
  }
  return(disp)
}

#' Phylogenetic/functional uniqueness of individuals or species.
#' @description Average dissimilarity between an individual or species and all others in a community with replacement.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree An hclust or phylo object.
#' @param abund A boolean (T/F) indicating whether dissimilarity should be calculated per individual (T) or species (F). If not specified, default is TRUE.
#' @param relative A boolean (T/F) indicating whether dissimilarity should be relative to the maximum distance between any two species in the tree.
#' @details Uniqueness is the originality measure of Pavoine et al. (2005).
#' @return A matrix of sites x species values.
#' @references Pavoine, S., Ollier, S. & Dufour, A.-B. (2005) Is the originality of a species measurable? Ecology Letters, 8: 579-586.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' uniqueness(comm, tree)
#' uniqueness(comm, tree, abund = FALSE)
#' uniqueness(comm, tree, abund = FALSE, relative = TRUE)
#' @export
uniqueness <- function(comm, tree, abund = TRUE, relative = FALSE){
  if(!abund)
    comm <- ifelse(comm > 0, 1, 0)
  if(missing(tree))
    tree = hclust(as.dist(matrix(1,ncol(comm),ncol(comm))))
  comm <- as.matrix(comm)
  unique <- matrix(0,nrow(comm),ncol(comm))

  for (i in 1:nrow(comm)){    					                         #cycle through all sites/samples
    present <- which(comm[i,]>0)                                 #which species exist in this site
    nSpp <- length(present)                                      #how many species are present in this site
    proportion <- comm[i,present]/sum(comm[i,])                  #proportion incidence/abundance of species in this site
    distance <- as.matrix(cophenetic(tree))[present,present]     #cophenetic distances of species in this site
    for (r in 1:nSpp)
      unique[i,present[r]] <- sum(sum(distance[r,] * proportion))
  }
  if(relative)
    unique <- unique / max(cophenetic(tree))
	return(unique)
}

#' Scaled mean squared error of accumulation curves.
#' @description Accuracy (scaled mean squared error) of accumulation curves compared with a known true diversity value (target).
#' @param accum A matrix resulting from the alpha.accum or beta.accum functions (sampling units x diversity values).
#' @param target The true known diversity value, with which the curve will be compared. If not specified, default is the diversity observed with all sampling units.
#' @details Among multiple measures of accuracy (Walther & Moore 2005) the SMSE presents several advantages, as it is (Cardoso et al. 2014):
#' (i) scaled to true diversity, so that similar absolute differences are weighted according to how much they represent of the real value;
#' (ii) scaled to the number of sampling units, so that values are independent of sample size;
#' (iii) squared, so that small, mostly meaningless fluctuations around the true value are down-weighted; and
#' (iv) independent of positive or negative deviation from the real value, as such differentiation is usually not necessary.
#' For alpha diversity accuracy may also be weighted according to how good the data is predicted to be. The weight of each point in the curve is proportional to its sampling intensity (i.e. n/Sobs).
#' @return Accuracy values (both raw and weighted) for all observed and estimated curves.
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, in press.
#' @references Walther, B.A. & Moore, J.L. (2005) The concepts of bias, precision and accuracy, and their use in testing the performance of species richness estimators, with a literature reviewof estimator performance. Ecography, 28, 815-829.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' acc.alpha = alpha.accum(comm1)
#' accuracy(acc.alpha)
#' accuracy(acc.alpha, 10)
#' acc.beta = beta.accum(comm1, comm2, tree)
#' accuracy(acc.beta)
#' accuracy(acc.beta, c(1,1,0))
#' @export
accuracy <- function(accum, target = -1){
  if(ncol(accum) > 5 || accum[nrow(accum), 3] > 1){		#if alpha
		if (target == -1)
			target <- accum[nrow(accum), 3]
		intensTotal = accum[nrow(accum), 2] / accum[nrow(accum), 3]	#sampling intensity = final n / final S
		if(ncol(accum) > 10){              #if non-parametric
      smse <- matrix(0, 13, nrow = 2)
      for (i in 1:nrow(accum)){
      	intensity = accum[i, 2] / accum[i, 3] / intensTotal
      	error = (accum[i,3] - target)^2 / (target^2 * nrow(accum))
      	smse[1,1] <- smse[1,1] + error
      	smse[2,1] <- smse[2,1] + error * intensity
      	for (j in 2:13){
      		error = (accum[i,j+6] - target)^2 / (target^2 * nrow(accum))
      		smse[1,j] <- smse[1,j] + error
      		smse[2,j] <- smse[2,j] + error * intensity
      	}
      }
      rownames(smse) <- c("Raw", "Weighted")
     	colnames(smse) <- c("Obs", "Jack1ab", "Jack1abP", "Jack1in", "Jack1inP", "Jack2ab", "Jack2abP", "Jack2in", "Jack2inP", "Chao1", "Chao1P", "Chao2", "Chao2P")
    }
    else{                              #if curve
      smse <- matrix(0, 5, nrow = 2)
      for (i in 3:nrow(accum)){
      	intensity = accum[i, 2] / accum[i, 3] / intensTotal
      	for (j in 1:5){
          if (!is.na(accum[i,j+2])){
          	error = (accum[i,j+2] - target)^2 / (target^2 * nrow(accum))
          	smse[1,j] <- smse[1,j] + error
          	smse[2,j] <- smse[2,j] + error * intensity
          }
      	}
      }
      rownames(smse) <- c("Raw", "Weighted")
      colnames(smse) <- c("Obs", "Clench", "Exponential", "Rational", "Weibull")
    }
	} else {																						#if beta
		if (target[1] == -1)
			target <- accum[nrow(accum), 2:4]
		smse <- rep(0, 3)
		for (i in 1:nrow(accum)){
			for (j in 1:3)
				smse[j] <- smse[j] + (accum[i,j+1] - target[j])^2
		}
		smse <- smse / nrow(accum)
		smse <- list(Btotal=smse[1], Brepl=smse[2], Brich=smse[3])
		smse <- c(unlist(smse))
	}
	return(smse)
}

#' Slope of accumulation curves.
#' @description This is similar to the first derivative of the curves at each of its points.
#' @param accum A matrix resulting from the alpha.accum or beta.accum functions (sampling units x diversity values).
#' @details Slope is the expected gain in diversity when sampling a new individual. The slope of an accumulation curve, of either observed or estimated diversity, allows verifying if the asymptote has been reached (Cardoso et al. 2011).
#' This is an indication of either the completeness of the inventory (low final slopes of the observed curve indicate high completeness) or reliability of the estimators (stability of the slope around a value of 0 along the curve indicates reliability).
#' @return A matrix of sampling units x slope values.
#' @references Cardoso, P., Pekar, S., Jocque, R. & Coddington, J.A. (2011) Global patterns of guild composition and functional diversity of spiders. PLoS One, 6, e21710.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' acc.alpha = alpha.accum(comm1)
#' slope(acc.alpha)
#' acc.beta = beta.accum(comm1, comm2, tree)
#' slope(acc.beta)
#' @export
slope <- function(accum){
	if(ncol(accum) > 5 || accum[nrow(accum), 3] > 1){			#if alpha
		sl <- accum[,-2]
		accum <- rbind(rep(0,ncol(accum)), accum)
		for (i in 1:nrow(sl)){
			sl[i,1] <- i
			for (j in 2:ncol(sl)){
				sl[i,j] <- (accum[i+1,j+1]-accum[i,j+1])/(accum[i+1,2]-accum[i,2])
			}
		}
	} else {																							#if beta
		sl <- accum
		sl[1,] <- 0
		sl[1,1] <- 1
		for (i in 2:nrow(sl)){
			for (j in 2:ncol(sl)){
				sl[i,j] <- (accum[i,j]-accum[i-1,j])
			}
		}
	}
	return(sl)
}

#' Optimization of alpha diversity sampling protocols.
#' @description Optimization of alpha diversity sampling protocols when different methods and multiple samples per method are available.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only to optimize PD or FD sampling).
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param base A vector defining a base protocol from which to build upon (complementarity analysis) (length must be equal to number of methods).
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @param prog Present a text progress bar in the R console.
#' @details Often a combination of methods allows sampling maximum plot diversity with minimum effort, as it allows sampling different sub-communities, contrary to using single methods.
#' Cardoso (2009) proposed a way to optimize the number of samples per method when the target is to maximize sampled alpha diversity. It is applied here for TD, PD and FD, and for one or multiple sites simultaneously.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix of samples x methods (values being optimum number of samples per method). The last column is the average alpha diversity value, rescaled to 0-1 if made for several sites, where 1 is the true diversity of each site.
#' @references Cardoso, P. (2009) Standardization and optimization of arthropod inventories - the case of Iberian spiders. Biodiversity and Conservation, 18, 3949-3962.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2), c(4,3,2))
#' colnames(comm) <- c("Sp1","Sp2","Sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' optim.alpha(comm,,methods)
#' optim.alpha(comm, tree, methods)
#' optim.alpha(comm,, methods = methods, base = c(0,0,1), runs = 100)
#' @export
optim.alpha <- function(comm, tree, methods, base, runs = 1000, prog = TRUE){
	
	##preliminary stats
	methods <- as.vector(t(methods))
	nSamples <- length(methods)							##number of samples
	metUnique <- as.vector(t(unique(methods)))				##list of methods
	metNum <- length(metUnique)							##number of methods
	if (missing(base))										##if no samples to start with for complementarity analysis
		samples <- rep(0,metNum)
	else
		samples <- base
	nMiss <- nSamples - sum(samples)				##number of samples missing
	nSamplesMet <- rep(0,metNum)						##samples per method
	for (m in 1:metNum)											
		nSamplesMet[m] <- sum(methods == metUnique[m])
	
	##accumulation process
	if (prog) pb <- txtProgressBar(max = nMiss+1, style = 3)
	div <- rep(0,nMiss+1)										##diversity along the optimal accumulation curve
	if (sum(samples) > 0)															
		div[1] <- optim.alpha.stats(comm, tree, methods, samples, runs)
	if (prog) setTxtProgressBar(pb, 1)
	for (s in 2:(nMiss+1)){
		samples <- rbind (samples, rep(0,metNum))
		samples[s,] <- samples[s-1,]
		metValue <- rep(0, metNum)										#diversity when adding each method
		for (m in 1:metNum){
			if (samples[s,m] < nSamplesMet[m]){
				samples[s,m] <- samples[s,m] + 1
				metValue[m] <- optim.alpha.stats(comm, tree, methods, samples[s,], runs)
				samples[s,m] <- samples[s,m] - 1
			}
		}
		div[s] <- max(metValue)
		best <- which(metValue == div[s])
		if (length(best) > 1)
			best = best[sample(1:length(best),1)]						#if tie, choose one of the best methods randomly
		samples[s, best] <- samples[s, best] + 1
		if (prog) setTxtProgressBar(pb, s)
	}
	if (prog) close(pb)
	colnames(samples) <- metUnique
	rownames(samples) <- (0:nMiss+sum(samples[1,]))
	samples <- cbind(samples, div)
	return(samples)
}

#' Efficiency statistics for alpha-sampling.
#' @description Average alpha diversity observed with a given number of samples per method.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only to optimize PD or FD sampling).
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param samples A vector defining the number of samples per method to be evaluated (length must be equal to number of methods).
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @details Different combinations of samples per method allow sampling different sub-communities.
#' This function allows knowing the average TD, PD or FD values for a given combination, for one or multiple sites simultaneously.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A single average alpha diversity value. Rescaled to 0-1 if made for several sites, where 1 is the true diversity of each site.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2), c(4,3,2))
#' colnames(comm) <- c("Sp1","Sp2","Sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' optim.alpha.stats(comm,,methods, c(1,1,1))
#' optim.alpha.stats(comm, tree, methods = methods, samples = c(0,0,1), runs = 100)
#' @export
optim.alpha.stats <- function(comm, tree, methods, samples, runs = 1000){
	
	##preliminary stats
	if (!missing(tree))
		tree <- xTree(tree)
	if(length(dim(comm)) == 3)					##number of sites
		nSites <- dim(comm)[3]
	else
		nSites <- 1
	methods <- as.vector(t(methods))
	metUnique <- as.vector(t(unique(methods)))				##list of methods
	metNum <- length(metUnique)					##number of methods
	div <- 0														##average diversity obtained using this particular combination of samples per method
	
	for (i in 1:nSites){
		if (nSites > 1){
			site <- as.matrix(comm[,,i])
			true <- sobs(site, tree) 				##true diversity of each site
		} else {
			site <- as.matrix(comm)
			true <- 1
		}
		
		for (r in 1:runs){
			addSample <- rep(0, ncol(comm))
			for (m in 1:metNum){
				if (samples[m] > 0){
					filterList <- site[which(methods == metUnique[m]),,drop=F]						##filter by method m
					filterList <- filterList[sample(nrow(filterList),samples[m]),,drop=F]	##randomly select rows
					addSample <- rbind(addSample, filterList)															##add random samples
				}
			}
			div <- div + sobs(addSample, tree) / runs / nSites / true
		}
	}
	return(div)
}

#' Optimization of beta diversity sampling protocols.
#' @description Optimization of beta diversity sampling protocols when different methods and multiple samples per method are available.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only to optimize PD or FD sampling).
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param base Allows defining a base mandatory protocol from which to build upon (complementarity analysis). It should be a vector with length = number of methods.
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @param prog Present a text progress bar in the R console.
#' @details Often, comparing differences between sites or the same site along time (i.e. measure beta diversity) it is not necessary to sample exhaustively. A minimum combination of samples targeting different sub-communities (that may behave differently) may be enough to perceive such differences, for example, for monitoring purposes.
#' Cardoso et al. (in prep.) introduce and differentiate the concepts of alpha-sampling and beta-sampling. While alpha-sampling optimization implies maximizing local diversity sampled (Cardoso 2009), beta-sampling optimization implies minimizing differences in beta diversity values between partially and completely sampled communities.
#' This function uses as beta diversity measures the Btotal, Brepl and Brich partitioning framework (Carvalho et al. 2012) and respective generalizations to PD and FD (Cardoso et al. 2014).
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix of samples x methods (values being optimum number of samples per method). The last column is the average absolute difference from real beta.
#' @references Cardoso, P. (2009) Standardization and optimization of arthropod inventories - the case of Iberian spiders. Biodiversity and Conservation, 18, 3949-3962.
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Cardoso, P., et al. (in prep.) Optimal inventorying and monitoring of taxon, phylogenetic and functional diversity.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm3 <- matrix(c(2,0,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2, comm3), c(4,3,3))
#' colnames(comm) <- c("sp1","sp2","sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' optim.beta(comm, methods = methods, runs = 100)
#' optim.beta(comm, tree, methods = methods, abund = TRUE, base = c(0,0,1), runs = 100)
#' @export
optim.beta <- function(comm, tree, methods, base, abund = FALSE, runs = 1000, prog = TRUE){
	
	##preliminary stats
	methods <- as.vector(t(methods))
	nSamples <- length(methods)							##number of samples
	metUnique <- as.vector(t(unique(methods)))				##list of methods
	metNum <- length(metUnique)							##number of methods
	
	if (missing(base))										##if no samples to start with
		samples <- rep(0,metNum)
	else
		samples <- base
	nMiss <- nSamples - sum(samples)				##number of samples missing
	nSamplesMet <- rep (0, metNum)					##samples per method
	for (m in 1:metNum)											
		nSamplesMet[m] <- sum(methods == metUnique[m])
	
	##accumulation process
	if (prog) pb <- txtProgressBar(max = nMiss+1, style = 3)
	diff <- rep(0,nMiss+1)														#absolute difference along the optimal accumulation curve
	diff[1] <- optim.beta.stats(comm, tree, methods, samples, abund, runs)
	if (prog) setTxtProgressBar(pb, 1)
	if (diff[1] == "NaN")
		diff[1] = 1
	for (s in 2:(nMiss+1)){
	  samples <- rbind (samples, rep(0,metNum))
		samples[s,] <- samples[s-1,]
		metValue <- rep(1, metNum)										#absolute difference when adding each method
		for (m in 1:metNum){
			if (samples[s,m] < nSamplesMet[m]){
				samples[s,m] <- samples[s,m] + 1
				metValue[m] <- optim.beta.stats(comm, tree, methods, samples[s,], abund, runs)
				samples[s,m] <- samples[s,m] - 1
			}
		}
		diff[s] <- min(metValue)
		best <- which(metValue == diff[s])
		if (length(best) > 1)
			best = best[sample(1:length(best),1)]						#if tie, choose one of the best methods randomly
		samples[s, best] <- samples[s, best] + 1
		if (prog) setTxtProgressBar(pb, s)
	}
	if (prog) close(pb)
	colnames(samples) <- metUnique
	rownames(samples) <- (0:nMiss+sum(samples[1,]))
	samples <- cbind(samples, diff)
	return(samples)
}

#' Efficiency statistics for beta-sampling.
#' @description Average absolute difference between sampled and real beta diversity when using a given number of samples per method.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree An hclust or phylo object (used only to optimize PD or FD sampling).
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param samples The combination of samples per method we want to test. It should be a vector with length = number of methods.
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @details Different combinations of samples per method allow sampling different sub-communities.
#' This function allows knowing the average absolute difference between sampled and real beta diversity for a given combination, for one or multiple sites simultaneously.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A single average absolute beta diversity difference value.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm3 <- matrix(c(2,0,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2, comm3), c(4,3,3))
#' colnames(comm) <- c("sp1","sp2","sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' optim.beta.stats(comm,,methods, c(1,1,1))
#' optim.beta.stats(comm, tree, methods = methods, samples = c(0,0,1), runs = 100)
#' @export
optim.beta.stats <- function(comm, tree, methods, samples, abund = FALSE, runs = 1000){

	##preliminary stats
	if(length(dim(comm)) == 3){					##number of sites
		nSites <- dim(comm)[3]
	}else{
		message("need sample data from at least two sites to perform analyses")
		return(0)
	}
	methods <- as.vector(t(methods))
	metUnique <- as.vector(t(unique(methods)))				##list of methods
	metNum <- length(metUnique)					##number of methods
	diff <- 0														##average absolute difference between observed and true diversity obtained using this particular combination of samples per method
	
	##calculate true beta values
	sumComm <- matrix(0, nrow = nSites, ncol = ncol(comm))
	for (i in 1:nSites){
		sumComm[i,] <- colSums(comm[,,i])
	}
	true <- beta(sumComm, tree, abund)
	
	##calculate absolute difference between sampled and true beta values
	for (r in 1:runs){
		sumComm <- matrix(0, nrow = nSites, ncol = ncol(comm))
		for (m in 1:metNum){
			if (samples[m] > 0){
				filterList <- comm[which(methods == metUnique[m]),,,drop=F] 							##filter by method m
				filterList <- filterList[sample(nrow(filterList),samples[m]),,,drop=F]		##randomly select rows
				for (i in 1:nSites){
					sumComm[i,] <- sumComm[i,] + colSums(filterList[,,i,drop=F])
				}
			}
		}
		sampleBeta <- beta(sumComm, tree, abund)
		for(i in 1:3){
			diff <- diff + mean(abs(sampleBeta[[i]] - true[[i]])) / 3 / runs
		}
	}
	return(diff)
}

#' Optimization of spatial sampling.
#' @description Optimization of sampling site distribution in space based on environmental (or other) variables.
#' @param layers A Raster* object (typically a multi-layer type: RasterStack or RasterBrick).
#' @param n The number of intended sampling sites (clusters).
#' @param latlong Boolean indicating whether latitude and longitude should be taken into account when clustering.
#' @param clusterMap Boolean indicating whether to build a new raster with clusters.
#' @details Optimizing the selection of sampling sites often requires maximizing the environmental diversity covered by them.
#' One possible solution to this problem, here adopted, is performing a k-means clustering using environmental data and choosing the sites closest to the multidimensional environmental centroid of each cluster for sampling (Jimenez-Valverde & Lobo 2004)
#' @return Either a matrix of cells x clusters (also indicating distance to centroid, longitude and latitude of each cell) or a list with such matrix plus the clusterMap.
#' @references Jimenez-Valverde, A., & Lobo, J. M. (2004). Un metodo sencillo para seleccionar puntos de muestreo con el objetivo de inventariar taxones hiperdiversos: el caso practico de las familias Araneidae y Thomisidae (Araneae) en la comunidad de Madrid, Espana. Ecologia, 18: 297-305.
#' @export
optim.spatial <- function(layers, n, latlong = TRUE, clusterMap = TRUE){
	dataMat <- as.matrix(layers)
	dataMat <- dataMat[complete.cases(dataMat),]
	dataMat <- cbind(dataMat, rasterToPoints(layers[[1]])[,1:2])					##add latlong
	if (latlong)
		res <- kmeans(dataMat, n)        ##do k-means
	else
		res <- kmeans(dataMat[,-c((ncol(dataMat)-1),ncol(dataMat))], n)        ##do k-means

	cl = c()
	for(c in 1:n){
		cData <- dataMat[res$cluster==c,]           #filter to cluster c
		cCenter <- res$centers[c,]
		dist2centroid <- c()
		for(r in 1:nrow(cData))
			dist2centroid[r] = dist(rbind(cData[r,], cCenter))
		cData <- cbind(rep(c, nrow(cData)), dist2centroid, cData[,ncol(cData)], cData[,(ncol(cData)-1)])
		colnames(cData) <- c("cluster", "dist2centroid", "lat", "long")
		cData <- cData[sort.list(cData[,2]), ]
		cl <- rbind(cl, cData)
	}

	#output raster with clusters
	if(clusterMap){
		map <- rasterize(rasterToPoints(layers[[1]])[,1:2], layers[[1]], res$cluster)
		names(map) <- "clusters"
		cl <- list(cl, map)
	}
	return(cl)
}

#' Species-area relationship (SAR).
#' @description Fits and compares several of the most supported models for the species (or PD, or FD) -area relationship.
#' @param comm Either a vector with the diversity values per site, or a sites x species matrix.
#' @param tree An hclust or phylo object (used only to fit the PD or FD-area relationships, requires comm to be a sites x species matrix).
#' @param area A vector with the area per site.
#' @details Larger areas (often islands) usually carry more species. Several formulas were proposed in the past to describe this relationship (Arrhenius 1920, 1921; Gleason 1922).
#' Recently, the same approach began to be used for other measures of diversity, namely phylogenetic (PD) and functional (FD) diversity (Whittaker et al. 2014).
#' The function compares some of the most commonly used and theoretically or empirically suported models.
#' The relationships for PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix with the different model parameters and explanatory power.
#' @references Arrhenius, O. (1920) Distribution of the species over the area. Meddelanden fran Vetenskapsakadmiens Nobelinstitut, 4: 1-6.
#' @references Arrhenius, O. (1921) Species and area. Journal of Ecology, 9: 95-99.
#' @references Gleason, H.A. (1922) On the relation between species and area. Ecology, 3: 158-162.
#' @references Whittaker, R.J., Rigal, F., Borges, P.A.V., Cardoso, P., Terzopoulou, S., Casanoves, F., Pla, L., Guilhaumon, F., Ladle, R. & Triantis, K.A. (2014) Functional biogeography of oceanic islands and the scaling of functional diversity in the Azores. Proceedings of the National Academy of Sciences USA, 111: 13709-13714.
#' @examples div <- c(1,2,3,4,4)
#' comm <- matrix(c(2,0,0,0,3,1,0,0,2,4,5,0,1,3,2,5,1,1,1,1), nrow = 5, ncol = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:4), method="euclidean"), method="average")
#' area <- c(10,40,80,160,160)
#' sar(div,,area)
#' sar(comm,,area)
#' sar(comm,tree,area)
#' @export
sar <- function(comm, tree, area){
	if(is.vector(comm)){
		div = comm
	} else if (missing(tree)){
		div = alpha(comm)
	} else {
		div = alpha(comm, tree)
	}
	results <- matrix(NA, 6, 7)
	colnames(results) <- c("c", "z", "r2", "AIC", "\U0394 AIC", "AICc", "\U0394 AICc")
	rownames(results) <- c("Linear", "Linear (origin)", "Exponential", "Exponential (origin)", "Power", "Power (origin)")
	k <- c(3,2,3,2,3,2)
	model <- list()
	model[[1]] <- try(nls(div ~ c + z*area, start = data.frame(c = 0, z = 1)))
	model[[2]] <- try(nls(div ~ z*area, start = data.frame(z = 1)))
	model[[3]] <- try(nls(div ~ c + z*log(area), start = data.frame(c = 0, z = 1)))
	model[[4]] <- try(nls(div ~ z*log(area), start = data.frame(z = 1)))
	model[[5]] <- try(nls(div ~ c + area^z, start = data.frame(c = 0, z = 1)))
	model[[6]] <- try(nls(div ~ area^z, start = data.frame(z = 1)))
	for(m in 1:length(model)){
		if(k[m] == 3){
			results[m,1] <- coef(summary(model[[m]]))[1,1]
			results[m,2] <- coef(summary(model[[m]]))[2,1]
		} else {
			results[m,2] <- coef(summary(model[[m]]))[1,1]
		}
		pred <- predict(model[[m]], area=area)
		results[m,3] <- r2(pred, div)
		results[m,4] <- AIC(pred, div, k[m])
		results[m,6] <- AICc(pred, div, k[m])
	}
	for(m in 1:length(model)){
		results[m,5] <- results[m,4] - min(results[,4])
		results[m,7] <- results[m,6] - min(results[,6])
	}
	return(results)
}	

#' General dynamic model of oceanic island biogeography (GDM).
#' @description Fits and compares several of the most supported models for the GDM (using TD, PD or FD).
#' @param comm Either a vector with the diversity values per island, or an island x species matrix.
#' @param tree An hclust or phylo object (used only to fit the PD or FD GDM, requires comm to be a sites x species matrix).
#' @param area A vector with the area of islands.
#' @param time A vector with the age of islands. If not given, the species-area relationship is returned instead.
#' @details The general dynamic model of oceanic island biogeography was proposed to account for diversity patterns within and across oceanic archipelagos as a function of area and age of the islands (Whittaker et al. 2008).
#' Several different equations have been found to describe the GDM, extending the different SAR models with the addition of a polynomial term using island age and its square (TT2), depicting the island ontogeny.
#' The first to be proposed was an extension of the exponential model (Whittaker et al. 2008), the power model extensions following shortly after (Fattorini 2009; Steinbauer et al. 2013), as was the linear model (Cardoso et al. subm.).
#' The relationships for PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix with the different model parameters and explanatory power.
#' @references Cardoso, P., Borges, P.A.V., Carvalho, J.C., Rigal, F., Gabriel, R., Cascalho, J. & Correia, L. (subm.) Automated discovery of relationships, models and principles in ecology. Pre-print available from bioRxiv doi: http://dx.doi.org/10.1101/027839
#' @references Fattorini, S. (2009) On the general dynamic model of oceanic island biogeography. Journal of Biogeography, 36: 1100-1110.
#' @references Steinbauer, M.J, Klara, D., Field, R., Reineking, B. & Beierkuhnlein, C. (2013) Re-evaluating the general dynamic theory of oceanic island biogeography. Frontiers of Biogeography, 5: 185-194.
#' @references Whittaker, R.J., Triantis, K.A. & Ladle, R.J. (2008) A general dynamic theory of oceanic island biogeography. Journal of Biogeography, 35: 977-994.
#' @examples div <- c(1,3,5,8,10)
#' comm <- matrix(c(2,0,0,0,3,1,0,0,2,4,5,0,1,3,2,5,1,1,1,1), nrow = 5, ncol = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:4), method="euclidean"), method="average")
#' area <- c(10,40,80,160,160)
#' time <- c(1,2,3,4,5)
#' gdm(div,,area,time)
#' gdm(comm,tree,area,time)
#' gdm(div,,area)
#' @export
gdm <- function(comm, tree, area, time){
	if(missing(time))
		return(sar(comm,tree,area))
	if(is.vector(comm)){
		div = comm
	} else if (missing(tree)){
		div = alpha(comm)
	} else {
		div = alpha(comm, tree)
	}
	results <- matrix(NA, 4, 9)
	colnames(results) <- c("c", "z", "x", "y", "r2", "AIC", "\U0394 AIC", "AICc", "\U0394 AICc")
	rownames(results) <- c("Linear", "Exponential", "Power (area)", "Power (area, time)")
	k <- 5
	model <- list()
	model[[1]] <- try(nls(div ~ c + z*area + x*time + y*time^2, start = data.frame(c=0, z=1, x=1, y=0)))
	model[[2]] <- try(nls(div ~ c + z*log(area) + x*time + y*time^2, start = data.frame(c=0, z=1, x=1, y=0)))
	model[[3]] <- try(nls(div ~ exp(c + z*log(area) + x*time + y*time^2), start = data.frame(c=0, z=1, x=1, y=0)))
	model[[4]] <- try(nls(div ~ exp(c + z*log(area) + x*log(time) + y*log(time)^2), start = data.frame(c=0, z=1, x=1, y=0)))
	for(m in 1:length(model)){
		results[m,1] <- coef(summary(model[[m]]))[1,1]
		results[m,2] <- coef(summary(model[[m]]))[2,1]
		results[m,3] <- coef(summary(model[[m]]))[3,1]
		results[m,4] <- coef(summary(model[[m]]))[4,1]
		pred <- predict(model[[m]], area=area, time=time)
		results[m,5] <- r2(pred, div)
		results[m,6] <- AIC(pred, div, k)
		results[m,8] <- AICc(pred, div, k)
	}
	for(m in 1:length(model)){
		results[m,7] <- results[m,6] - min(results[,6])
		results[m,9] <- results[m,8] - min(results[,8])
	}
	return(results)
}	

#' Interspecific abundance-occupancy relationship (IAOR).
#' @description Fits and compares several of the most supported models for the IAOR.
#' @param comm A sites x species matrix with abundance values.
#' @details Locally abundant species tend to be widespread while locally rare species tend to be narrowly distributed.
#' That is, for a given species assemblage, there is a positive interspecific abundance-occupancy relationship (Brown 1984).
#' This function compares some of the most commonly used and theoretically or empirically suported models (Nachman 1981; He & Gaston 2000; Cardoso et al. subm.).
#' @return A matrix with the different model parameters and explanatory power.
#' @references Brown, J.H. (1984) On the relationship between abundance and distribution of species. American Naturalist, 124: 255-279.
#' @references Cardoso, P., Borges, P.A.V., Carvalho, J.C., Rigal, F., Gabriel, R., Cascalho, J. & Correia, L. (subm.) Automated discovery of relationships, models and principles in ecology. Pre-print available from bioRxiv doi: http://dx.doi.org/10.1101/027839
#' @references He, F.L. & Gaston, K.J. (2000) Estimating species abundance from occurrence. American Naturalist, 156: 553-559.
#' @references Nachman, G. (1981) A mathematical model of the functional relationship between density and spatial distribution of a population. Journal of Animal Ecology, 50: 453-460.
#' @examples comm <- matrix(c(4,3,2,1,5,4,3,2,3,2,1,0,6,3,0,0,0,0,0,0), nrow = 5, ncol = 4, byrow = TRUE)
#' iaor(comm)
#' @export
iaor <- function(comm){
  results <- matrix(NA, 4, 7)
  colnames(results) <- c("a", "b", "r2", "AIC", "\U0394 AIC", "AICc", "\U0394 AICc")
  rownames(results) <- c("Linear", "Exponential", "Negative Binomial", "SR")
  k <- c(3,3,2,2)
  abund <- colMeans(comm)                   #mean abundance per species (including sites with 0 individuals)
  occup <- colMeans(ifelse(comm>0,1,0))     #proportion occupancy per species

  model <- list()
  model[[1]] <- try(nls(logit(occup) ~ a+b*log(abund), start = data.frame(a = 1, b = 1))) #linear
  model[[2]] <- try(nls(occup ~ 1-exp(a*abund^b), start = data.frame(a = -1, b = 1))) #exponential
  model[[3]] <- try(nls(occup ~ 1-(1+(abund/a))^(0-a), start = data.frame(a = 0))) #negative binomial
  model[[4]] <- try(nls(occup ~ abund/(a+abund), start = data.frame(a = 0))) #SR = Clench with asymptote 1
  for(m in 1:length(model)){
    if(m < 3){
      results[m,1] <- coef(summary(model[[m]]))[1,1]
      results[m,2] <- coef(summary(model[[m]]))[2,1]
    } else {
      results[m,1] <- coef(summary(model[[m]]))[1,1]
    }
    pred <- predict(model[[m]], abund=abund)
    if(m==1) pred = revLogit(pred)
    results[m,3] <- r2(pred, occup)
    results[m,4] <- AIC(pred, occup, k[m])
    results[m,6] <- AICc(pred, occup, k[m])
  }
  for(m in 1:length(model)){
    results[m,5] <- results[m,4] - min(results[,4])
    results[m,7] <- results[m,6] - min(results[,6])
  }
  return(results)
}

#' Simulation of species abundance distributions (SAD).
#' @description Creates artificial communities following given SADs.
#' @param n total number of individuals.
#' @param s number of species.
#' @param sad The SAD distribution type (lognormal, uniform, broken stick or geometric). Default is lognormal.
#' @param sd The standard deviation of lognormal distributions. Default is 1.
#' @details Species Abundance Distributions may take a number of forms. A lognormal SAD probably is the most supported by empirical data, but we include other common types useful for testing multiple algorithms including several of the functions in BAT. 
#' @return A matrix of species x abundance per species.
#' @examples comm1 <- sim.sad(10000, 100)
#' comm2 <- sim.sad(10000, 100, sd = 2)
#' comm3 <- sim.sad(10000, 100, sad = "uniform")
#' par(mfrow=c(1,3))
#' hist(log(comm1$Freq))
#' hist(log(comm2$Freq))
#' hist(log(comm3$Freq))
#' @export
sim.sad <- function(n, s, sad = "lognormal", sd = 1) {
	if (s > n)
		stop("Number of species can't be larger than number of individuals")
	sppnames = paste("Sp", 1:s, sep="") ##species names
	sad <- match.arg(sad, c("lognormal", "uniform", "broken", "geometric"))
	
	##lognormal distribution
	switch(sad, lognormal = {
		comm = sample(sppnames, size = n, replace = T, prob = c(rlnorm(s, sdlog = sd)))
		
		##uniform distribution
	}, uniform = {
		comm = sample(sppnames, size = n, replace = T)
		
		##broken stick distribution
	}, broken = {
		broken.stick <- function(p){
			result = NULL
			for(j in 1:p) {
				E = 0
				for(x in j:p)
					E = E+(1/x)
				result[j] = E/p
			}
			return(result)
		}
		broken.prob = broken.stick(s)
		comm = sample(sppnames, size = n, replace = TRUE, prob = c(broken.prob))
	}, geometric = {
		geo.ser <- function(s, k = 0.3){
			result = NULL
			for (x in 1:s) {
				result[x] = k*(1-k)^(x-1)/(1-(1-k)^s)
			}
			return(result)
		}
		geo.prob = geo.ser(s)
		comm = sample(sppnames, size = n, replace = TRUE, prob = c(geo.prob))
	})
	return(as.data.frame(table(comm)))
}

#' Simulation of species spatial distributions.
#' @description Creates artificial communities with given SAD and spatial clustering.
#' @param n total number of individuals.
#' @param s number of species.
#' @param sad The SAD distribution type (lognormal, uniform, broken stick or geometric). Default is lognormal.
#' @param sd The standard deviation of lognormal distributions. Default is 1.
#' @param dist The spatial distribution of individual species populations (aggregated, random, uniform or gradient). Default is aggregated.
#' @param clust The clustering parameter (higher values create more clustered populations). Default is 1.
#' @details The spatial distribution of individuals of given species may take a number of forms.
#' Competitive exclusion may cause overdispersion, specific habitat needs or cooperation may cause aggregation and environmental gradients may cause abundance gradients.
#' @return A matrix of individuals x (species, x coords and y coords).
#' @examples par(mfrow = c(3 ,3))
#' comm = sim.spatial(100, 9, dist = "uniform")
#' for(i in 1:9){
#' 	sp <- comm[comm[1] == paste("Sp", i, sep = ""), ]
#' 	plot(sp$x, sp$y, main = paste("Sp", i), xlim = c(0,1), ylim = c(0,1))
#' }
#' comm = sim.spatial(1000, 9, sad = "lognormal", sd = 0.5, dist = "aggregated", clust = 2)
#' for(i in 1:9){
#' 	sp <- comm[comm[1] == paste("Sp", i, sep=""), ]
#' 	plot(sp$x, sp$y, main = paste("Sp", i), xlim = c(0,1), ylim = c(0,1))
#' }
#' @export
sim.spatial <- function(n, s, sad = "lognormal", sd = 1, dist = "aggregated", clust = 1){
	repeat{
		simsad <- sim.sad(n, s, sad, sd)
		dist <- match.arg(dist, c("aggregated", "random", "uniform", "gradient"))
	
		##aggregated distribution
		switch(dist, aggregated = {
			clust <- 1/(4*clust)
			cluster <- vector("list", s)
			for (i in 1:s)
				cluster[[i]] = rThomas(1, clust, n, as.owin(c(0,1,0,1)))
			ppx <- NULL
			ppy <- NULL
			for (j in 1:s){
				ppx <- c(ppx, cluster[[j]]$x[1:simsad[,2][j]])
				ppy <- c(ppy, cluster[[j]]$y[1:simsad[,2][j]])
			}
			spp <- rep(as.character(simsad[,1]), simsad[,2])
			comm <- data.frame(Spp = spp, x = ppx, y = ppy)
			
			##random distribution
		}, random = {
			rand <- runifpoint(n, as.owin(c(0,1,0,1)))
			spp <- rep(as.character(simsad[,1]), simsad[,2])
			comm <- data.frame(Spp = spp, x = rand$x, y = rand$y)
			
			##uniform distribution
		}, uniform = {
			rand <- rSSI(1/n, n, as.owin(c(0,1,0,1)))
			spp <- rep(as.character(simsad[,1]), simsad[,2])
			comm <- data.frame(Spp = spp, x = rand$x, y = rand$y)
			
			##gradient distribution
		}, gradient = {
			comm <- rpoint(n, function(x,y){x})
			spp <- rep(as.character(simsad[,1]), simsad[,2])
			comm <- data.frame(Spp = spp, x = comm$x, y = comm$y)
		})
		if(nrow(comm) == n && length(which(is.na(comm$x))) == 0){
			break
		}
	}
	return(comm)
}

#' Plots of simulated species spatial distributions.
#' @description Plots individuals from artificial communities with given SAD and spatial clustering.
#' @param comm artificial community data from function sim.spatial.
#' @param sad boolean indicating if the SAD plot should also be shown. Default is FALSE.
#' @param s number of species to plot simultaneously. Default is the number of species in comm.
#' @details Function useful for visualizing the results of sim.spatial.
#' @examples comm <- sim.spatial(1000, 24)
#' sim.plot(comm)
#' sim.plot(comm, sad = TRUE)
#' sim.plot(comm, s = 9)
#' @export
sim.plot <- function(comm, sad = FALSE, s = 0){
	spp <- length(unique(comm$Spp))
	if(s < 1){
		if(!sad)
			side = ceiling(sqrt(spp))
		else
			side = ceiling(sqrt(spp+1))
	}	else{
		side = ceiling(sqrt(s))
	}
	par(mfrow = c(side,side), mar = c(1,1,2,1))
	if(sad)
		hist(log(table(comm[,1])), main = "All species", xlab = "Abundance (log)", xaxt = "n")
	for(i in 1:spp){
		sp <- comm[comm[1] == paste("Sp", i, sep=""), ]
		plot(sp$x, sp$y, main = paste("Sp", i), xlim = c(0,1), ylim = c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")
	}
}

#' Simulation of sampling from artificial communities.
#' @description Simulates a sampling process from artificial communities.
#' @param comm simulated community data from function sim.spatial.
#' @param cells number of cells to divide the simulated space into. Default is 100.
#' @param samples number of samples (cells) to randomly extract. Default is the number of cells (the entire community).
#' @details The space will be divided in both dimensions by sqrt(cells).
#' @details Function useful for simulating sampling processes from the results of sim.spatial.
#' @details May be used as direct input to other functions (e.g. alpha, alpha.accum, beta, beta.accum) to test the behavior of multiple descriptors and estimators.
#' @return A matrix of samples x species (values are abundance per species per sample).
#' @examples comm <- sim.spatial(1000, 10)
#' sim.sample(comm)
#' sim.sample(comm, cells = 10, samples = 5)
#' @export
sim.sample <- function(comm, cells = 100, samples = 0){
	side <- round(sqrt(cells),0)
	cells = side^2

	comm$ind <- 0
	xv <- cut(comm$x, seq(0, 1, 1/side))
	yv <- cut(comm$y, seq(0, 1, 1/side))
	grid1 <- data.frame(table(xv, yv))
	grid1 <- grid1[,-3]
	
	s <- 1:cells
	for (i in 1:cells){
		id <- NULL
		id <- which (xv == grid1$xv[s[i]] & yv == grid1$yv[s[i]])
		comm$ind[id] <- paste("Sample", i, sep="")
	}
	comm <- table(comm$ind, comm$Spp) ##entire community
	comm <- comm[rownames(comm) != "0", ]
	if (samples < 1 || samples > nrow(comm))
		samples = nrow(comm)
	
	##number of samples to take
	samp <- comm[sample(nrow(comm), samples, replace = FALSE),] ## sampled community
	
	return(samp)
}

#' Simulation of phylogenetic or functional tree.
#' @description Simulates a random tree.
#' @param s number of species.
#' @param m a structural parameter defining the average difference between species. Default is 100. Lower numbers create trees dominated by increasingly similar species, higher numbers by increasingly dissimilar species.
#' @details A very simple tree based on random genes/traits.
#' @return An hclust object.
#' @examples tree <- sim.tree(10)
#' plot(as.dendrogram(tree))
#' tree <- sim.tree(100,10)
#' plot(as.dendrogram(tree))
#' tree <- sim.tree(100,1000)
#' plot(as.dendrogram(tree))
#' @export
sim.tree <- function(s, m = 100){
	sim.matrix <- matrix(sample(0:m, ceiling(s*m/50), replace = TRUE), nrow = s, ncol = m)
	tree <- hclust(dist(sim.matrix), method = "average")
	tree$height <- tree$height / max(tree$height)
	return(tree)
}

#' Sample data of spiders in Arrabida (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 sampling units. Details are described in:
#' Cardoso, P., Gaspar, C., Pereira, L.C., Silva, I., Henriques, S.S., Silva, R.R. & Sousa, P. (2008) Assessing spider species richness and composition in Mediterranean cork oak forests. Acta Oecologica, 33: 114-127.
#' 
#' @docType data
#' @keywords datasets
#' @name arrabida
#' @usage data(arrabida)
#' @format A data frame with 320 sampling units (rows) and 338 species (variables).
NULL

#' Sample data of spiders in Geres (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 sampling units. Details are described in:
#' Cardoso, P., Scharff, N., Gaspar, C., Henriques, S.S., Carvalho, R., Castro, P.H., Schmidt, J.B., Silva, I., Szuts, T., Castro, A. & Crespo, L.C. (2008) Rapid biodiversity assessment of spiders (Araneae) using semi-quantitative sampling: a case study in a Mediterranean forest. Insect Conservation and Diversity, 1: 71-84.
#' 
#' @docType data
#' @keywords datasets
#' @name geres
#' @usage data(geres)
#' @format A data frame with 320 sampling untis (rows) and 338 species (variables).
NULL

#' Sample data of spiders in Guadiana (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 sampling units. Details are described in:
#' Cardoso, P., Henriques, S.S., Gaspar, C., Crespo, L.C., Carvalho, R., Schmidt, J.B., Sousa, P. & Szuts, T. (2009) Species richness and composition assessment of spiders in a Mediterranean scrubland. Journal of Insect Conservation, 13: 45-55.
#' 
#' @docType data
#' @keywords datasets
#' @name guadiana
#' @usage data(guadiana)
#' @format A data frame with 192 sampling units (rows) and 338 species (variables).
NULL

#' Functional tree for 338 species of spiders
#'
#' A dataset representing the functional tree for 338 species of spiders captured in Portugal.
#' For each species were recorded: average size, type of web, type of hunting, stenophagy, vertical stratification in vegetation and circadial activity. Details are described in:
#' Cardoso, P., Pekar, S., Jocque, R. & Coddington, J.A. (2011) Global patterns of guild composition and functional diversity of spiders. PLoS One, 6: e21710.
#' 
#' @docType data
#' @keywords datasets
#' @name functree
#' @usage data(functree)
#' @format An hclust object with 338 species.
NULL

#' Taxonomic tree for 338 species of spiders (surrogate for phylogeny)
#'
#' A dataset representing an approximation to the phylogenetic tree for 338 species of spiders captured in Portugal.
#' The tree is based on the linnean hierarchy, with different suborders separated by 1 unit, families by 0.75, genera by 0.5 and species by 0.25.
#' 
#' @docType data
#' @keywords datasets
#' @name phylotree
#' @usage data(phylotree)
#' @format An hclust object with 338 species.
NULL