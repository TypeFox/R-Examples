calcmpprob <- function(object, chr, step, mapfx, ibd, mrkpos)
{
	# if step < 0, values are computed at the midpoints of intervals only
	# if step = 0, values are computed only at markers
	if (mapfx=="kosambi") mf <- kosambiX2R  else mf <- haldaneX2R

	n.founders <- nrow(object$founders)
	n.finals <- nrow(object$finals)

	#map for the probabilities
	nmap <- list()
	prob <- list()
	for (ii in chr)
	{
		chrmap <- object$map[[ii]] 
		lmap <- length(chrmap)

		#vector of indicies, which convert markers for this chromosome into columns of the data matrix
		mrkchr <- match(names(chrmap), colnames(object$finals))
		#discard indicies relating to markers which have no data
		mrkchr <- mrkchr[!is.na(mrkchr)]

		# set up the grid positions at which to compute probs
		if (step<0)
		{
			#if step < 0 then we're only computing probabilities at midpoints, so it's pretty simple
			r12 <- r23 <- sapply(diff(chrmap), mf)
			pos <- chrmap[1:(length(chrmap)-1)]+diff(chrmap)/2
			
			#left marker of this point
			left <- 1:length(r12)
			#right marker of this point
			right <- 2:length(chrmap)
			###modification, unnecessary line deleted
			mid <- rep(-1, length(r12))
		}
		if (step==0)
	 	{
		pos <- chrmap
		## this will take care of all middle markers
		r12 <- r23 <- mf(diff(chrmap))
		r12 <- c(r12[1], r12[1:(length(r12)-1)], r12[length(r12)-1])
		r23 <- c(r23[2], r23[2:length(r23)], r23[length(r23)])
		## need to do endpoints separately. 
		left <- c(1, 1:(length(chrmap)-2), length(chrmap)-2)
		mid <- left+1
		right <- mid+1
		}
		if (step>0) 
		{
			###modification, this should start at min(chrmap), not 0, if we start at zero there may not be flanking markers.
			pos <- seq(min(chrmap), max(chrmap), step)
			## if mrkpos==TRUE, include map positions as well
			if (mrkpos)  pos <- sort(unique(c(pos, chrmap)))
			#which of these poisitions are markers?
			pinc <- pos %in% chrmap
			#for the ones which are markers, preserve the names
			names(pos)[pinc] <- names(chrmap)[which(chrmap %in% pos)]
			#otherwise names them loc <position in cM>
			names(pos)[!(pinc)] <- paste("loc", pos[!pinc], sep="")
			#convert the positions which are markers to indicies of markers
			flag <- match(pos, chrmap)
			#for those that aren't markers, set the flag to 0
			flag[is.na(flag)] <- 0
			
			#cumsum(flag > 0) gives a vector where the element i counts the number of markers in positions 1:i. That is, the first value will be
			#a 1 because the first position is a marker.
			###modification. This was calculating the probability at marker 2 as the first probability, and the probability at the second last marker last, regardless of where those should actually have fit in the output list of probabilities.
			#cumsum(flag > 0) gives the number of markers to the left of the current position (excluding the current position, if it is a marker). Add 1 because it will start at 0 (the first value of flag is 0)
			left <- cumsum(flag > 0)[1:(length(flag)-2)]
			#add on (flag > 0) because we want to step the right marker across one more if the current position is a marker
			right <- left + (flag > 0)[2:(length(flag)-1)] + 1
			r12 <- sapply(((pos[2:(length(pos)-1)]-chrmap[left])), mf)
			r23 <- sapply(((chrmap[right]-pos[2:(length(pos)-1)])), mf)	
			mid = flag[2:(length(flag)-1)]
			
			#add on the first and last markers, the left right and mid values are a bit different for these as they're treated differently in the
			#C code
			left = c(1, left, length(chrmap)-2)
			mid = c(2,mid, length(chrmap) - 1)
			right = c(3, right, length(chrmap))
			r12 = c(
				mf((chrmap[2] - chrmap[1])),r12,
				mf((chrmap[length(chrmap)-1] - chrmap[length(chrmap)-2]))
				)
			r23 = c(
				mf((chrmap[3] - chrmap[2])),r23,
				mf((chrmap[length(chrmap)] - chrmap[length(chrmap)-1]))
				)
		}
		if (mapfx=="kosambi")	r13 <- (r12+r23)/(1+4*r12*r23) else
		r13 <- r12+r23-2*r12*r23
		nmap[[ii]] <- pos

		n.pos <- length(r12)

		if (ibd==TRUE) {
		finals <- cbind(object$id, object$ibd$finals[,mrkchr])
			founders <- object$ibd$founders[,mrkchr]
		} else {
		finals <- cbind(object$id, object$finals[,mrkchr])
			founders <- object$founders[,mrkchr]
		}

		finals[is.na(finals)] <- -1

		finvec <- vector()
		fouvec <- vector()
		finals[,1] <- as.numeric(as.character(finals[,1]))

		for (i in 1:n.finals) 	 
		  finvec <- c(finvec, as.numeric(finals[i,]))

		for (i in 1:n.founders) fouvec <- c(fouvec, as.numeric(founders[i,]))  

		hap3pt <- .C("hap3ptfull",  as.integer(finvec), as.integer(fouvec), as.integer(object$pedigree[,1]), as.integer(object$pedigree[,2]), as.integer(object$pedigree[,3]), as.integer(n.finals), as.integer(n.founders), as.integer(lmap), as.integer(nrow(object$pedigree)), out=double(length=n.pos*n.founders*n.finals), as.double(r12), as.double(r23), as.double(r13), as.integer(left), as.integer(mid), as.integer(right), as.integer(n.pos), PACKAGE="mpMap")

		haps <- matrix(data=hap3pt$out, nrow=n.finals, ncol=n.founders*n.pos, byrow=TRUE) 
		colnames(haps) = paste("Position ", rep(pos, each=n.founders), ", Founder ", 1:n.founders, sep="")
		prob[[ii]] <- haps
	}
  attr(prob, "map") <- nmap 
  return(prob)
}
