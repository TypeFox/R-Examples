generate_obs <- function(geno, map, full.prob, fg, transpos, transval, founderld)
{
  # geno is formatted as a matrix. $finals, $founders, etc. 
  # will need to overlay on both founders and finals 
  obsgeno <- geno
  n.founders <- nrow(geno$founders)

  if (is.null(fg)) {
	  ## generate which markers will be biallelic
    o <- sample(c(TRUE, FALSE), ncol(geno$founders)/2, replace=TRUE, 
   	prob=c(full.prob, 1-full.prob))

    #	Could be an option to change later - allow for different founder prob
    #	strains
    #	  if (missing(domprob)) 
    domprob <- c(rep(1/(n.founders-1), n.founders-1), 0)
    dom <- matrix(data=0, nrow=n.founders, ncol=ncol(geno$founders)/2)

    ## founderld == TRUE means founder genotypes generated wrt rf
    if (founderld) {
	  ## create genotypes for founders at first biallelic markers
	  ndom <- sample(1:n.founders, size=1, prob=domprob)
	  dom[sample(1:n.founders, ndom), which(o==FALSE)[1]] <- 1

	  ### now create a genetic map for the markers where o==FALSE and 
	  ### adjacent recombination fractions
	  allmrk <- unlist(lapply(map, names))
	  dommrk <- allmrk[which(o==FALSE)]
	  chr <- rep(1:length(map), unlist(lapply(map, length)))
	  domchr <- chr[which(o==FALSE)]
	  pos <- unlist(map)
	  dompos <- pos[which(o==FALSE)]
	
	  recfr <- haldaneX2R(abs(diff(dompos)))
	  recfr[which(diff(domchr)!=0)] <- 0.5

	  ## now at each position generate a uniform according to recfr
	  for (i in 1:n.founders) {
		runi <- runif(length(recfr))
		rec <- (runi < recfr)
		rec <- c(dom[i, which(o==FALSE)[1]], rec)
		dom[i, which(o==FALSE)] <- cumsum(rec) %% 2
	  }
    	  o <- c(o,o)
	  dom <- cbind(dom, dom)

	  obsgeno$founders[,o==FALSE] <- dom[,o==FALSE]
    } else {
	# note, probably want to allow a general distribution for this
	ndom <- sample(1:n.founders, sum(o==FALSE), replace=TRUE, prob=domprob)
	ndom <- c(ndom, ndom)
	for (i in 1:sum(o==FALSE)) 
	  dom[sample(1:n.founders, ndom[i]),which(o==FALSE)[i]] <- 1
	dom <- cbind(dom, dom)
	
	obsgeno$founders[,o==FALSE] <- dom[,o==FALSE]
 
	### set special founders for translocation region if markers are biallelic
	if (transval>0) {
	  for (i in transpos) { 
	  if (o[i]==FALSE)
	#    obsgeno$founders[,i] <- rep(1,n.founders)
	     obsgeno$founders[transval,i] <- 0
	    }
	}
    } 

  for (i in 1:n.founders)
    obsgeno$finals[,obsgeno$founders[i,]==0][geno$finals[,obsgeno$founders[i,]==0]==i] <- 0
  obsgeno$finals[,o==FALSE][(obsgeno$finals[,o==FALSE]>0)] <- 1

  } else { ## else if founder genotypes are input
    ## since we're not taking out the qtl first, need to modify fg
    ## so that QTL founder genotypes are included
    allmrk <- unlist(lapply(map, names))
    chr <- rep(1:length(map), unlist(lapply(map, length)))
    pos <- unlist(map)

    recfr <- haldaneX2R(abs(diff(pos)))
    recfr[which(diff(chr)!=0)] <- 0.5

    qtlpos <- grep("QTL", allmrk)

    fg2 <- matrix(nrow=nrow(fg), ncol=ncol(fg)+length(qtlpos))
    colnames(fg2) <- allmrk
    rownames(fg2) <- rownames(fg)
    fg2[, match(colnames(fg), colnames(fg2))] <- as.matrix(fg)

    ## what do you do if you get a false? have to change the input genotypes
    ## which doesn't make a lot of sense. have to assume for this setup that
    ## QTL genotype matches whichever flanking recfr is smaller. 
    matchpos <- qtlpos
    for (i in 1:length(qtlpos)) 
	matchpos[i] <- (qtlpos[i]+1)*(recfr[qtlpos[i]]<recfr[qtlpos[i]-1]) + (qtlpos[i]-1)*(recfr[qtlpos[i]]>=recfr[qtlpos[i]-1])
    fg2[, qtlpos] <- fg2[,matchpos]

    obsgeno$founders <- cbind(fg2, fg2)
    for (i in 1:ncol(obsgeno$finals))
	obsgeno$finals[,i] <- obsgeno$founders[,i][obsgeno$finals[,i]]
  }
  return(obsgeno)
}
# end function: generate.obs

