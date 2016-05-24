partitionCountData <-
function(data, n.boots.part=100, max.pop=500){
	
# convert data format to table of individuals
	ind.data <- classToInd(data)
	
	c <- ind.data[which(ind.data[,2] != 0), 2]
	c <- sort(c, decreasing=FALSE)
	
	if (length(unique(c)) == 1){
		return(print("All capture counts are the same. Partitioning not possilbe."))
		
	}
	
	
	# --- internal function ---
	get.partition.vectors <- function(cut.1, cut.2, n.inds, set){
		part.1 <- NA; part.2 <- NA; part.3 <- NA
		if (cut.2>0){
			part.1 <- seq(1,cut.2)
			}
		if (cut.2==0){part.1 <- NA}
	
		if (cut.1<n.inds){
			part.3 <- seq(cut.1, n.inds)
			}
		if (cut.1==n.inds){
			if (cut.2==0){
				part.3 <- seq(1,n.inds)
				} else{
				part.3 <- NA
				}
			}
		part.1.3 <- union(part.1, part.3)
		nas <- which(is.na(part.1.3))
		if (length(nas)>0){part.1.3 <- part.1.3[-nas]}
		part.2 <- seq(1, n.inds)
		part.2 <- part.2[-part.1.3]
		if (length(part.2)==0){part.2 <- NA}
	
		count.1 <- set[part.1]
		count.2 <- set[part.2]
		count.3 <- set[part.3]
		parts <- list(count.1, count.2, count.3, part.1, part.2, part.3)
		names(parts) <- c("count.1", "count.2", "count.3", "part.1", "part.2", "part.3")
		return(parts)
		}
	
	# --- internal function ---
	find.3.way.partitions <- function(set){
		n.inds <- length(set)
		poss.parts <- matrix(nrow=n.inds^2, ncol=2)
		part.i <- 0
		for (i in n.inds:1){
			part.1 <- n.inds
			remaining <- seq(1, i,by=1)
			for (j in max(remaining):0){
				part.i <- part.i + 1
				poss.parts[part.i,] <- c(i,j) 
				}
			}
		n.parts <- min(which(is.na(poss.parts[,1]))) - 1
		poss.parts <- poss.parts[1:n.parts,]
		ln.L.each <- rep(NA, n.parts)
		count.part.v <- rep(NA, n.parts)
		for (part.i in 1:n.parts){
			cut.1 <- poss.parts[part.i, 1]
			cut.2 <- poss.parts[part.i, 2]
			results <-get.partition.vectors(cut.1, cut.2, n.inds, set)
			count.1 <- results$count.1
			count.2 <- results$count.2
			count.3 <- results$count.3
			part.1 <- results$part.1
			part.2 <- results$part.2
			part.3 <- results$part.3
		
			bin.1 <- is.na(part.1[1])==FALSE
			bin.2 <- is.na(part.2[1])==FALSE
			bin.3 <- is.na(part.3[1])==FALSE
			count.part.v[part.i] <- sum(bin.1, bin.2, bin.3)
		 
			p.1 <- 1/length(part.1)
			p.2 <- 1/length(part.2)
			p.3 <- 1/length(part.3)
			if (is.na(part.1[1])==FALSE){
				prob.1 <- dmultinom(count.1, prob=rep(p.1, length(part.1)), log=TRUE)
			} else{prob.1 <- NA}
				if (is.na(part.2[1])==FALSE){
					prob.2 <- dmultinom(count.2, prob=rep(p.2, length(part.2)), log=TRUE)
			} else{prob.2 <- NA}
				if (is.na(part.3[1])==FALSE){
					prob.3 <- dmultinom(count.3, prob=rep(p.3, length(part.3)), log=TRUE)
			} else{prob.3 <- NA}
		
			ln.L.each[part.i] <- sum(prob.1, prob.2, prob.3,  na.rm=TRUE)
			}
	
		summary.m <- cbind(poss.parts, ln.L.each, count.part.v)
		colnames(summary.m) <- c("cut.1", "cut.2", "ln.L", "n.groups")
		summary.m <- as.data.frame(summary.m)
		return(summary.m)
	}
	
	# --- internal function ---
	get.possible.partitions <- function(set){
		n.inds <- length(set)
		poss.parts <- matrix(nrow=n.inds^2, ncol=2)
		part.i <- 0
		for (i in n.inds:1){
			part.1 <- n.inds
			remaining <- seq(1, i,by=1)
			for (j in max(remaining):0){
				part.i <- part.i + 1
				poss.parts[part.i,] <- c(i,j) 
				}
			}
		n.parts <- min(which(is.na(poss.parts[,1]))) - 1
		poss.parts <- poss.parts[1:n.parts,]
		return(poss.parts)	
		}
		
	# --- internal function ---
	delta.L.for.three.way.partitioning <- function(three.way.h){
		one.group <- which(three.way.h$n.groups==1)
		two.groups <- which(three.way.h$n.groups==2)
		three.groups <-  which(three.way.h$n.groups==3)
		ln.L.unparted <- three.way.h$ln.L[1]
		max.L.two.groups.ID <- two.groups[which.max(three.way.h$ln.L[two.groups])]
		max.L.two.groups.LIKE <-three.way.h$ln.L[max.L.two.groups.ID]
		max.L.three.groups.ID <- three.groups[which.max(three.way.h$ln.L[three.groups])]
		max.L.three.groups.LIKE <-three.way.h$ln.L[max.L.three.groups.ID]
		d.L.1.2 <- max.L.two.groups.LIKE - ln.L.unparted
		d.L.2.3 <-  max.L.three.groups.LIKE -  max.L.two.groups.LIKE
		d.L.1.3 <-  max.L.three.groups.LIKE -  ln.L.unparted
		delta.L <- list(d.L.1.2, d.L.2.3, d.L.1.3)
		names(delta.L) <- c("d.L.1.2", "d.L.2.3", "d.L.1.3") 
		return(delta.L)
		}
			
	# === Here we actually do control the analysis ===
	
	n.inds <- length(c)
	three.way <- find.3.way.partitions(c)
	
	one.group <- which(three.way$n.groups==1)
	two.groups <- which(three.way$n.groups==2)
	three.groups <-  which(three.way$n.groups==3)
	ln.L.unparted <- three.way$ln.L[1]
	max.L.two.groups.ID <- two.groups[which.max(three.way$ln.L[two.groups])]
	max.L.two.groups.LIKE <-three.way$ln.L[max.L.two.groups.ID]
	max.L.three.groups.ID <- three.groups[which.max(three.way$ln.L[three.groups])]
	max.L.three.groups.LIKE <-three.way$ln.L[max.L.three.groups.ID]
	
	p.2 <- max.L.two.groups.ID
	poss.parts <- get.possible.partitions(c)
	cut.1 <- poss.parts[p.2, 1]
	cut.2 <- poss.parts[p.2, 2]
	results <-get.partition.vectors(cut.1, cut.2, n.inds, c)
	bin.1 <- as.numeric(is.na(results$part.1[1])==FALSE)
	bin.2 <- as.numeric(is.na(results$part.2[1])==FALSE)
	bin.3 <- as.numeric(is.na(results$part.3[1])==FALSE)
	
	bins <- c(bin.1, bin.2, bin.3)
	two.parts <- which(bins==1)	
	counts <- list(results$count.1, results$count.2, results$count.3)
	count.2.A <- unlist(counts[two.parts[1]])
	count.2.B <- unlist(counts[two.parts[2]])
	
	p.3 <- max.L.three.groups.ID
	poss.parts <- get.possible.partitions(c)
	cut.1 <- poss.parts[p.3, 1]
	cut.2 <- poss.parts[p.3, 2]
	results <-get.partition.vectors(cut.1, cut.2, n.inds, c)
	count.3.A <- results$count.1
	count.3.B <- results$count.2	
	count.3.C <- results$count.3
	low.2.thirds <- c(count.3.A, count.3.B)
	up.1.third <- count.3.C

	# --- These are the observed differences in log-likelihoods
	d.L.1.2 <- max.L.two.groups.LIKE - ln.L.unparted
	d.L.2.3 <-  max.L.three.groups.LIKE -  max.L.two.groups.LIKE 
	d.L.1.3 <-  max.L.three.groups.LIKE - ln.L.unparted
	
	fitT <- suppressWarnings(fitTirm(data, max.pop))
	
	boot.d.L.2.3.TIRM <- rep(NA, n.boots.part)
	for (boot.i in 1:n.boots.part){
		sim.data <- simTirm(na=fitT$ml.na, nb=fitT$ml.nb, alpha=fitT$alpha, s=fitT$sample.size)
		boot.set.AB <- classToInd(sim.data)[,2]
		boot.results <- find.3.way.partitions(boot.set.AB)
		boot.d.L.2.3.TIRM[boot.i] <- delta.L.for.three.way.partitioning(boot.results)$d.L.2.3
	}
	p.2.3.TIRM <- length(which(boot.d.L.2.3.TIRM >= d.L.2.3)) / (n.boots.part + 1)
	p.2.3 <- p.2.3.TIRM
	
	low.2.thirds <- ind.data[-which(ind.data[,2] %in% up.1.third), c(1,2)]
	
	up.1.third <- ind.data[which(ind.data[,2] %in% up.1.third),c(1,2)]
		
## convert outputs back to class format
	
	low.2.thirds <- indToClass(low.2.thirds)
	
	up.1.third <- indToClass(up.1.third)
		
	
	results <- list(low.2.thirds, up.1.third, p.2.3)
	names(results) <- c("low.2.thirds", "up.1.third", "p.2.3")
	return(results)
	
}
