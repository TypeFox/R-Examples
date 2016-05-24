## wrapper function for calculating matrix based resemblance measures
## in neighborhoods, that means multi-plot similarity or dissimilarity 
## indices that are calculated for each plot in a (regular) array and 
## its neighbors. Which surrounding plots make up the neighborhood
## can be defined with giving either a radius or two radiuses (that than
## define a ring).
"rin" <-
function(veg, coord, dn, func, test=TRUE, permutations=100, permute=2, sfno=TRUE, p.level=0.05, ...){
	p <- permutations
	p.l <- p.level
	veg <- data.frame(veg)
	coord <- data.frame(coord)
	if (any(rownames(veg)!=rownames(coord))){
		warning("Plot names did not conform between species matrix and coordinates!")
	}
	## function that's later needed internally in the case of permute = 3
	# permutation is that species in the focal plot are randomly selected
	# p times from the pool (with sfno=TRUE, species from the neighborhood only,
	# with sfno=FALSE, species from veg)
	# the function works plot wise
	if(sfno){
	momo <- function(dat){
			vek.rndm <- lapply(c(1:p), function(x) with(dat, as.numeric(sample(x[foc,colSums(x)>0]))))
			foc <- with(dat, which(rownames(x)==foc))
			vek.rndm <- lapply(vek.rndm, function(y) rbind(y, dat$x[-foc,colSums(dat$x)>0]))
			sim.perms <- sapply(vek.rndm, function(x) eval(parse(text=func)))
			return(sim.perms)
		}
	}
	else{
	momo <- function(dat){
			vek.rndm <- lapply(c(1:p), function(x) with(dat, as.numeric(sample(x[foc,]))))
			foc <- with(dat, which(rownames(x)==foc))
			vek.rndm <- lapply(vek.rndm, function(y) rbind(y, dat$x[-foc,]))
			sim.perms <- sapply(vek.rndm, function(x) eval(parse(text=func)))
			return(sim.perms)
		}	
	}
	
	## function to calculate matrix based resemblance measures in neighborhoods
	"rin.in" <- function(veg, coord, dn, func, permute=2, ...){
	rownames(veg) <- c(1:nrow(veg))
    rownames(coord) <- rownames(veg)
    plots <- rownames(veg)
    dist.all <- dist(coord[,1:2])
    if(is.character(dn)){
    	dist.nbs <- apply(as.matrix(dist.all), 2, rank, ties.method="random")
    	dist.nbs <- ifelse(dist.nbs <= as.numeric(dn), 1, 0)
    }
    else{
    	if(max(dist.all) < min(dn)){
        	stop("Are you sure that the neighbor definition is correct?")
    	}
    	if(length(dn)==1){
        	dist.nbs <- ifelse(as.matrix(dist.all) <= dn, 1, 0)
        	dimnames(dist.nbs) <- list(plots, plots)
    	}
    	else{
        	dist.nbs <- ifelse(((as.matrix(dist.all) == 0) | ((as.matrix(dist.all) >= min(dn)) & (as.matrix(dist.all) <= max(dn)))), 1, 0)
        	dimnames(dist.nbs) <- list(plots, plots)
    	}
    }
    # which number of plots in a neighborhood is valid (depends on resemblance coefficient)
    # methods that treat the focal plot special need more plots in the neighborhood
	what.valid <- grepl("foc=foc", func) + 1
	#calculate the number of plots in each neighborhood
   	n.plots <- rowSums(dist.nbs)
	# make a selector vector for the valid neighborhood data sets
	sub.valid <- n.plots > what.valid
	#correct the entrys in the n.plots vector
   	n.plots <- n.plots[sub.valid]
	# create a list of neighborhood data sets
	tmp <- lapply(plots[sub.valid], function(x) list(foc=x, x=veg[dist.nbs[,x]==1,]))
   	# calculate the number of species in each neighborhood
   	n.spec <- sapply(tmp, function(x) sum(colSums(x[[2]])>0))
   	
   	
   	## now it diverts. 
   	# for permute = 1 | 2 mps per neighborhood is calculated for whole matrix
   	if(permute < 3) {
   		dis <- sapply(tmp, function(x) with(x, eval(parse(text=func)))) 
    	dis.mat <- data.frame(n.plots = rep(NA, length(plots)), n.spec = rep(NA, length(plots)), dis = rep(NA, length(plots)))
    	dis.mat[sub.valid,] <- cbind(n.plots, n.spec, dis)
    	names(dis.mat)[3] <- names(dis)[1]
		res <- dis.mat
	}	
	# for permute = 3
	else {
##	gegebenenfalls hier noch tmp aufbereiten, auf dass nur noch die arten
##	als pool berücksichtigt werden die in der jeweiligen neighborhood vorkommen!
		sim.perm.pre <- t(sapply(tmp, function(x) momo(x)))
		sim.perms <- matrix(nrow=length(plots), ncol=p)
    	sim.perms[sub.valid,] <- sim.perm.pre
		res <- sim.perms
	}
	return(res)
	}

	## calculate matrix based resemblance measure in neighborhoods
	sim.start <- rin.in(veg, coord=coord, dn=dn, func = func, ...)
	if(test){
		## permute the original species matrix (permute across rows)
		if (permute == 2){
			veg.rndm <- lapply(c(1:p), function(x) apply(veg, permute, function(x) sample(x)))
			# calculate do mps for each permuted matrix, write results into list
			sim.perms <- sapply(veg.rndm, function(x) rin.in(x, coord=coord, dn=dn, func=func, ...)[,3])
		}
		## permute the original species matrix (permute across cols)
		else if (permute == 1){
			veg.rndm <- lapply(c(1:p), function(x) t(apply(veg, permute, function(x) sample(x))))
			# calculate do mps for each permuted matrix, write results into list
			sim.perms <- sapply(veg.rndm, function(x) rin.in(x, coord=coord, dn=dn, func=func, ...)[,3])	
		}
		## permute the species in the focal plot (drawn randomly from the
		## available species pool)
		else if (permute == 3){
			sim.perms <- rin.in(veg, coord=coord, dn=dn, func=func, permute=3, ...)
		}
		
		## do the significance stuff
		sig.mat <- sim.perms - sim.start[,3]
		sig.high <- (rowSums(sig.mat > 0)+1)/(p+(1/p))
		sig.low <- (rowSums(sig.mat < 0)+1)/-(p+(1/p))
		## compare against the rowMeans of the results of the permuted runs, 
		## if it is above mean, significance is tested on the upper tail 
		## and if it is below mean, significance is tested on the lower tail
		sim.mean <- apply(sim.perms, 1, "mean")
		sim.test <- ifelse(sim.start[,3] < sim.mean, sig.low, sig.high)
		sig <- ifelse(abs(sim.test) < p.l, "*", "ns")
		sig.prefix <- ifelse(sim.test < 0, "-", "+")
		sims <- cbind(sim.start, data.frame(p.val = sim.test, sig = sig, sig.sign = sig.prefix))
		## achtung, nur für testzwecke.
		## sims <- list(sims, sim.perms)
	}
	else{
		sims <- sim.start
	}
	return(sims)
}
