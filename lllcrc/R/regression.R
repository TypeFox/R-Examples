#' Local averaging for LLLMs
#' 
#' A nearest-neighbors procedure is used in conjunction with the Epanechnikov
#' kernel to define a kernel smooth of multinomial outcomes across the
#' covariate space
#' 
#' See Kurtz 2013, Chapter on multiple sclerosis
#' 
#' @param dat The capture-recapture data in the form that is returned by
#' \code{\link{formatdata}} or \code{\link{micro.post.stratify}}.
#' @param kfrac The approximate fraction of the data that is included in the
#' support of the kernel for the local averages.
#' @param bw A matrix a single column, with rownames that match the covariate
#' names in \code{dat}.  The values in the column are scalars that are used in
#' constructing distances between covariate vectors.  Raw differences are
#' divided by the corresponding scalars before being squared in the context of
#' a Euclidean metric.
#' @return A list containing the original data (\code{dat}), the smoothed data
#' (\code{hpi}), and the effective sample sizes (\code{ess}) for each local
#' average, or row, in the smoothed data
#' @author Zach Kurtz
#' @references Kurtz 2013
#' @export smooth.patterns
#' @import data.table
smooth.patterns = function(dat, kfrac, bw)
{
	J = .I = NULL #The ONLY purpose of this line is to avoid a note about no variable bindings for J and .I in R CMD check
		      #For the inspiration, see http://stackoverflow.com/a/8096882/2232265
	#require(data.table)
	#if(!is.matrix(dat)) stop("dat needs to be a matrix") # replaced below:
	if(!is.matrix(dat)) dat = data.matrix(dat)
	nd = nrow(dat)
	if(!is.element("mct", colnames(dat))) dat[,"mct"] = rep(1,nd)
	nc = sum(dat[,"mct"])
	x.covs = colnames(dat)[substr(colnames(dat), 1,1) == "x"]
	n.covs = length(x.covs)
	if(is.null(bw)) {
		bw = matrix(1, nrow = n.covs, ncol = 1)
	}else if(!is.matrix(bw)){
		bw = as.matrix(bw)
	}
	if(length(bw) != n.covs) stop("The length of the vector bw needs to equal the number of covariates x.[dis/con]...")
	y.bits = colnames(dat)[substr(colnames(dat), 1,1) == "y"]
	ny = length(y.bits)
	# Initialize hat_Pi, the smoothed capture-pattern frequencies,
	#	plus a column for effective sample sizes:
	hpi = matrix(0, ncol = ny, nrow = nd)
	colnames(hpi) = substring(y.bits, 2) #ess = effective sample size
	ess = matrix(0, ncol = 1, nrow = nd); colnames(ess) = "ess"
	covs = dat[,x.covs, drop = FALSE]

	for(i in 1:nd){
		dm = matrix(0, nrow = nd, ncol = n.covs) # Get the vector of distances
		for(j in 1:n.covs) dm[,j] = ((covs[,j] - covs[i,j])/bw[j])^2
		dvec = sqrt(rowSums(dm))
		# Identify the kernel boundary
		edf = cbind(dvec, dat[,"mct"], 1:nd) #empirical distribution of distances
		colnames(edf) = c("dvec", "mct", "origindex")
		edf = edf[order(edf[,1]),]
		ecdf = cumsum(edf[,"mct"])/nc #plot(1:nd, ecdf)
		dtt = data.table(ecdf, val = ecdf)
		setattr(dtt, "sorted", "ecdf") 
		edf.ind = dtt[J(kfrac), .I, roll = "nearest"][1]
		# (edf.ind is the highest index of a micro-postratum that lies within the kernel support)
		# Define weights only for points in the support (sup) of the kernel
		sup = edf[1:edf.ind,,drop = FALSE]
		if(sup[edf.ind,"dvec"] > 0) sup[,"dvec"] = sup[,"dvec"]/sup[edf.ind,"dvec"]
		wghts = (1.01-sup[,"dvec"]^2) # basic unscaled Epanechnikov weights
		wghts = wghts/sum(wghts*sup[,"mct"]) # normalize with replicates
		# Apply weights to get smooth capture pattern frequencies
		hpi[i,] = wghts %*% dat[sup[,"origindex", drop=FALSE],y.bits,drop=FALSE] #matrix(wghts, nrow = 1) 
		ess[i,1] = 1/sum(sup[,"mct"]*wghts^2 )
	}

	sdat = list(dat = dat, hpi = hpi, ess = ess)
	return(sdat)
}
