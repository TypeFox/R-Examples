####################################################################################
#' spotNormDesign
#'
#' Produces a normalized design and calculates the minimal distance 
#' if required.
#' A design  is a matrix with \code{dim} columns and \code{size} rows.
#' Distance can be calculated with respect to a fixed, nested design.
#' 
#' @param dim number, dimension of the problem (will be no. of columns of the result matrix)
#' @param size number of points with that dimension needed. (will be no. of rows of the result matrix).
#' @param calcMinDistance Boolean to indicate whether a minimal distance should be calculated.
#' @param nested nested design to be considered during distance calculation.
#' @param ineq_constr inequality constraint, smaller zero for infeasible points. Used to replace infeasible points with random points. Has to evaluate points in interval [0;1].
#'
#' @seealso This function is used as a basis for \code{\link{spotCreateDesignLhd}}.
#'
#' @return list \code{L}  \cr
#' - \code{L} consists of a matrix and nd (if required) a minimal distance 
#' @export
#' @keywords internal
#' @author Christian Lasarczyk, adapted for nested designs by Martin Zaefferer
####################################################################################
spotNormDesign <- function(dim,size, calcMinDistance=FALSE, nested=NULL, ineq_constr=NULL){
	step <- 1/size
	design <- replicate(dim, sample(0:(size-1),size) * step + runif(size) * step)

	if(!is.null(ineq_constr)){
		feasible <- apply(design,1,ineq_constr) <= 0
		design <- design[feasible,]
		while(nrow(design)<size){
			newP <- runif(dim)
			if(ineq_constr(newP)<=0){
				design<-rbind(design,as.numeric(newP))
			}
		}
	}
	
	des <- rbind(design, nested) #adds nested design points, to be considered in distance calculation.
	
	if (calcMinDistance)
		minDistance <- min(dist(des))
	else
		minDistance <- NA;
	
	list(design=design, minDistance=minDistance)
}


####################################################################################
#' spotCreateDesignLhd
#'
#' Create a latin hyper cube design based on the number of dimensions and the number of design points.
#' Designs are created repeatedly. The randomly created design with the lowest between-points-distance is chosen.
#' A given set of points can be considered during distance calculation (\code{spotConfig$nested.design}). This set may represent a nested, higher level design, as used for
#' optimization with Co-Kriging or similar approaches. Will be ignored if \code{spotConfig$nested.design} is NULL.
#'
#' @param spotConfig list of spotConfiguration
#' @param noDesPoints number of design points, default is NaN
#' @param retries number of retries, which is the number of trials to find a design with the lowest distance, default is 1
#'
#' @return matrix \code{design} \cr
#' - \code{design} has \code{dimension} columns and \code{noDesPoints} rows
#' with entries corresponding to the region of interest.
#' @export
#' @author Christian Lasarczyk, adapted for nested designs by Martin Zaefferer
#' @seealso \code{\link{spotCreateDesignBasicDoe}}, \code{\link{spotCreateDesignFrF2}}, 
#' \code{\link{spotCreateDesignLhs}}, \code{\link{spotCreateDesignLhsOpt}}
####################################################################################
spotCreateDesignLhd <- function(spotConfig, noDesPoints = NaN, retries= 1) {
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotCreateDesignLhd");
	#noDesPoints <- spotConfig$init.design.size
	## retries is the number of trials to find a 
	## design with the greates minimal distance, (default is 1)
	## Calculate required points for initial design
	
	## use roi or aroi:
	if(spotConfig$spot.fileMode){
		if(file.exists(spotConfig$io.aroiFileName))
			roiConfig <- spotReadRoi(spotConfig$io.aroiFileName,spotConfig$io.columnSep,spotConfig$io.verbosity)
		else
			roiConfig <- spotReadRoi(spotConfig$io.roiFileName,spotConfig$io.columnSep,spotConfig$io.verbosity)
	}else{
		roiConfig <- spotConfig$alg.aroi
		if(is.null(roiConfig)) roiConfig <- spotConfig$alg.roi
	}
	
	if(is.na(noDesPoints)) {
		dim <- nrow(roiConfig)
		## in last term, the second +1 is necessary for crossvalidation
		noDesPoints <- max(11*dim,1 + 3 * dim + dim * (dim - 1) / 2 + 1)
		## dim higher than 16, max will take the second arguent, else the first		
	}
	
	if(is.null(spotConfig$alg.constraint.ineq)){
		ineq_constr_norm <- NULL
	}else{
		ineq_constr_norm <- function(x){
			x <- roiConfig$lower + x * (roiConfig$upper - roiConfig$lower)
			spotConfig$alg.constraint.ineq(x)
		}
	}	
	
	## nested design: 
	nestdes <- spotConfig$nested.design
	## scale nested design to [0,1]
	for (param in row.names(roiConfig)){
		lowerBound <-  spotConfig$alg.roi[param,"lower"]
		upperBound <-  spotConfig$alg.roi[param,"upper"]
		nestdes[param] <- lowerBound + nestdes[param] * (upperBound-lowerBound)
	}
	
	## Bei einer Wiederholung muss die Distanz nicht berechnet werden
	best <- spotNormDesign(nrow(roiConfig),noDesPoints,calcMinDistance=retries>1,nested=nestdes,ineq_constr=ineq_constr_norm)
	
	if (retries>1) {
		for (i in 1:(retries-1)) {
			tmpDes <- spotNormDesign(nrow(roiConfig),noDesPoints,calcMinDistance=TRUE,nested=nestdes,ineq_constr=ineq_constr_norm)
			## maximize minimal distance
			if (tmpDes$minDistance > best$minDistance)
				best <- tmpDes
		}
	}
	
	design <- as.data.frame(best$design)
	colnames(design) <- row.names(roiConfig)
	for (param in row.names(roiConfig)){
		lowerBound <-  spotConfig$alg.roi[param,"lower"]
		upperBound <-  spotConfig$alg.roi[param,"upper"]
		
		## x.5 rounds to next even value: 2.25 -> 2.2 and 2.35 -> 2.4
		if (roiConfig[param,"type"]  == "INT" || roiConfig[param,"type"]  == "FACTOR"){
			## print( c(param, spotConfig$alg.roi[param,"type"], lowerBound))
			lowerBound <- lowerBound - 0.5
			upperBound <- upperBound + 0.4999999999999
		}
		design[param] <- lowerBound + design[param] * (upperBound-lowerBound)
	}
	
	## round Integers (and Factor)
	###if (any(spotConfig$alg.roi[["type"]] == "INT") || any(spotConfig$alg.roi[param,"type"]  == "FACTOR"))
	design[roiConfig[["type"]]  == "INT"] <- floor(design[roiConfig[["type"]]  == "INT"]+0.5)
	design[roiConfig[["type"]]  == "FACTOR"] <- floor(design[roiConfig[["type"]]  == "FACTOR"]+0.5)
	
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotCreateDesignLhd")
	design
}
