#' Compute founder probabilities for multi-parent crosses
#'
#' Using haplotype probabilities, computes the probability that each location on a genome was inherited from each founder. Locations are run either at markers only, at the midpoints of all intervals or at step sizes of x cM. Probabilities can be computed using internally, or with R/happy.hbrem or R/qtl.
#' @export mpprob
#' @aliases mpprob calcmpprob print.mpprob
#' @param object Object of class \code{mpcross}
#' @param chr Subset of chromosomes
#' @param step Step size (in cM) to create grid of positions at which to compute probabilities. At default value of 0, probabilities are calculated at marker positions only
#' @param threshold Threshold for calling founder probabilities
#' @param mapfx Map function used to convert map to recombination fractions
#' @param ibd Flag to indicate whether to compute probabilities using IBD genotypes
#' @param mrkpos Flag to indicate whether to compute probabilities at both marker positions and step size or just step size. Is overridden for step size of 0
#' @param program R package to use to compute probabilities
#' @param tempfiledirectory Directory in which to output temporary files. Default is current working directory
#' @param generations Number of generations to assume in HAPPY. see \code{happy}
#' @details If \code{program=="mpMap"} then probabilities are computed using flanking markers at positions across the genome and represent 3-point haplotype probabilities. If \code{program=="happy"} then probabilities are computed using default values in R/happy.hbrem, which calculates ancestral haplotypes without using pedigree information. This only allows for probabilities to be computed at midpoints of intervals. If \code{program=="qtl"} then probabilities are computed from multipoint founder probabilities in R/qtl. 
#' 
#' If \code{step<0} for R/mpMap or R/qtl or \code{step==0} for R/happy.hbrem, then probabilities are computed at the midpoints of marker intervals. However, if \code{step==0} for R/qtl or R/mpMap, probabilities are computed only at marker locations.
#' @return The input mpcross object is returned with two additional components:
#' \item{prob}{A list with founder probabilities for each chromosome. Format is a matrix with n.founders * n.markers columns and n.lines rows. Each group of n.founders columns will add up to 1. Founder probabilities are in the order of founders in the input founder matrix. }
#' \item{estfnd}{A list with estimated founders for each chromosome. Format is a matrix with n.markers columns and n.lines rows. Missing values indicate where no founder probability exceeded the input threshold. Numeric values for founders indicate the row in the input founder matrix corresponding to the estimated founder.}
#' @seealso \code{\link[mpMap]{plot.mpprob}}, \code{\link[mpMap]{summary.mpprob}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 50, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' plot(mpp.dat)
#' summary(mpp.dat)

mpprob <- function(object, chr, step=0, mrkpos=TRUE, mapfx=c("haldane", "kosambi"), ibd=FALSE, threshold=0.7, program=c("mpMap", "qtl", "happy"), tempfiledirectory="", generations=5)
{
  if (missing(mapfx)) mapfx="haldane"

  if (missing(chr))
  {
	chr <- names(object$map)
	print("No chromosomes specified, will default to all")
  }

  if (is.numeric(chr)) chr <- names(object$map)[chr]
  if (step<0) mrkpos <- FALSE ## override for midpoints with mpMap
  if (step==0) mrkpos <- TRUE ## override mrkpos flag in this case

  n.founders <- nrow(object$founders)
  n.finals <- nrow(object$finals)
  out <- object
 
  out$prob <- list()
  out$estfnd <- list()
  attr(out$estfnd, "threshold") <- threshold
  ## as a default, convert the probabilities to estimated founders as well.
  ffx <- function(x) 
	if (is.na(max(x)) | max(x) <= threshold) return(NA) else return(which(x>threshold))
 
  chr1 <- vector()
  for (i in chr) 
	if (length(object$map[[i]])==1) chr1 <- c(chr1, i)
  chr <- setdiff(chr, chr1)

  if (length(chr1)>0) 
	mrkpos <- TRUE

  if (program=="mpMap")
    out$prob <- calcmpprob(object, chr, step, mapfx, ibd, mrkpos)
  if (program=="qtl"){
	if (tempfiledirectory != "") 
	write2cross(object, filestem=paste(tempfiledirectory, "/tmp", sep=""))
	else write2cross(object, filestem="tmp")
   	cr <- qtl:::readMWril(tempfiledirectory, "tmp.ril.csv", "tmp.founder.csv", type=attr(object, "type"))
 	cr <- subset(cr, chr=chr)
	if (step >= 0)
	gp <- calc.genoprob(cr, step=step)
	else if (step < 0) {
		pos <- list()
		for (i in chr)
		{ 
		  pos[[i]] <- vector(length=length(object$map[[i]])*2-1)
		  pos[[i]][seq(1, length(pos[[i]]), 2)] <- object$map[[i]]
		  pos[[i]][seq(2, length(pos[[i]]), 2)] <- object$map[[i]][1:(length(object$map[[i]])-1)]+diff(object$map[[i]])/2
		  names(pos[[i]])[seq(1, length(pos[[i]]), 2)] <- names(object$map[[i]])
		  names(pos[[i]])[seq(2, length(pos[[i]]), 2)] <- paste("loc", round(pos[[i]][seq(2, length(pos[[i]]), 2)], 1), sep="")
		} 
		class(pos) <- "map"
	gp <- calc.genoprob2(cr, pos=pos)
	}
	if (length(chr1)>0)
	for (i in chr1) {
	m <- match(names(object$map[[i]]), colnames(object$finals))
	gp$geno[[i]]$prob <- array(dim=c(nrow(object$finals), 1, nrow(object$founders)))
	gp$geno[[i]]$prob[,1,] <- t(sapply(object$finals[,m], function(x) return(1*(x==object$founders[,m])/(sum(object$founders[,m]==x)))))
	attr(gp$geno[[i]]$prob, "map") <- object$map[[i]]
	}

	prob1 <- lapply(gp$geno, function(x) return(x$prob))
	prob <- lapply(prob1, function(x) {
		mat <- matrix(nrow=dim(x)[1], ncol=dim(x)[2]*dim(x)[3])
		for (i in 1:dim(x)[3])
			mat[, seq(i, ncol(mat), dim(x)[3])] <- as.matrix(x[,,i])
		return(mat)})
	crmap <- lapply(gp$geno, function(x) return(attr(x$prob, "map")))

	if (mrkpos==FALSE) {
		for (i in 1:length(prob)){
			if (step>0)
			m <- which(!(crmap[[i]] %in% seq(min(object$map[[i]]), max(object$map[[i]]), step)))
			else m <- which(names(crmap[[i]]) %in% names(object$map[[i]]))
			if (length(m)>0){
			m2 <- (rep(m, each=n.founders)-1)*n.founders+rep(1:n.founders, length(m))
			crmap[[i]] <- crmap[[i]][-m]
			prob[[i]] <- prob[[i]][,-m2]
			}
			colnames(prob[[i]]) <- paste("C",i,"P ", rep(1:length(crmap[[i]]), each=n.founders), ", Founder ", 1:n.founders, sep="")
			rownames(prob[[i]]) <- rownames(object$finals)
		}
	}

	attr(prob, "map") <- crmap

	out$prob <- prob
  }

  if (program=="happy"){
	require(happy.hbrem)
	prob <- list()
	map <- list()
	write2happy(object, filestem="tmp")
 	hin <- happy("tmp.data", "tmp.alleles", generations=generations, haploid=TRUE)
	# number of marker intervals
	nint <- unlist(lapply(object$map, length))-1
	cnint <- c(0, cumsum(nint))
	mrkint <- 1
	for (i in 1:length(chr)) {
	  prob[[i]] <- matrix(nrow=nrow(object$finals), ncol=nint[i]*n.founders)
	  while(mrkint <= cnint[i+1]) {
	    prob[[i]][, (mrkint-cnint[i]-1)*n.founders+1:n.founders] <- hdesign(hin, mrkint)
	    mrkint <- mrkint+1
	  }
  	  map[[i]] <- object$map[[i]][1:(length(object$map[[i]])-1)] + diff(object$map[[i]])/2
	  colnames(prob[[i]]) <- paste("C",i,"P ", rep(1:nint[i], each=n.founders), ", Founder ", 1:n.founders, sep="")
	  rownames(prob[[i]]) <- rownames(object$finals)
	} # end of chr loop
    class(map) <- "map"
    attr(prob, "map") <- map
    out$prob <- prob
  } # end of happy loop

  for (ii in c(chr, chr1))
  {
    haps <- out$prob[[ii]]    
    haps[is.nan(haps)] <- NA

    fmat <- matrix(nrow=nrow(haps), ncol=ncol(haps)/n.founders)
    for (kk in seq(1, ncol(haps), n.founders))
    {
 	fmat[,(kk-1)/n.founders+1] <- apply(haps[,kk:(kk+n.founders-1)], 1, ffx)
	fmat[,(kk-1)/n.founders+1] <- factor(fmat[,(kk-1)/n.founders+1], levels=1:n.founders)
    }
    out$estfnd[[ii]] <- fmat
  }

  names(out$prob) <- names(attr(out$prob,"map")) <- 
		names(out$estfnd) <- c(chr, chr1)
  
  attr(out$prob, "step") <- step
  attr(out$prob, "program") <- program
  attr(out$prob, "mapfx") <- mapfx
  attr(out$prob, "mrkpos") <- mrkpos

	#name all the rows of probabilities  
	for(i in names(out$prob))
	 	rownames(out$prob[[i]]) = paste("L", 1:n.finals, sep="")
	
  class(out) <- c("mpprob", class(object))
  return(out)
}
