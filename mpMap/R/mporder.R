#' Order markers within linkage groups
#' 
#' Orders markers within linkage groups using two-point or multipoint probabilities. Two-point ordering is based on estimated recombination fractions; multi-point ordering is based on R/qtl ripple function. 
#' @export
#' @param object Object of class \code{mpcross}
#' @param chr Selected chromosomes or linkage groups to order
#' @param type Which type of ordering to perform - two-point or multipoint
#' @param rm.rf Flag for whether to remove recombination fraction values from 2-pt ordering which have missing values
#' @param window Window size for multipoint ordering
#' @param repeats Number of times to repeat multipoint ordering
#' @param mapfx Map function to use to compute final cM positions
#' @param criterion Criterion used in 2-pt ordering to determine best order
#' @param missfx Function to use to fill missing recombination fractions. See \code{\link[mpMap]{fill}}
#' @param ... Additional arguments
#' @details \emph{Two-point ordering}\cr
#' To use the two-point ordering, the recombination fractions between all pairs of markers must first be estimated. If there are missing values in this matrix, the markers with the largest number of missing values will be removed until there are no missing values left. These markers will not be used in the ordering and are recommended to be inserted into the resulting framework map using \code{\link[mpMap]{add3pt}} later. 
#'
#' Multiple methods are used to investigate optimal two-point orderings. These are taken from the package \code{seriation} and include simulated annealing, hierarchical clustering, and traveling salesman solver. The orders are compared on the basis of the argument \code{criterion}. Thus the total path length, or sum of adjacent recombination fractions can be minimized; or the number of Anti-Robinson events/deviations; or the number of crossovers; or the sum of the adjacent two-point LOD scores. 
#'
#' \emph{Multi-point ordering} \cr
#' The multi-point ordering assumes that there is a pre-existing map, and then repeatedly applies the ripple function in R/qtl to investigate local permutations of the order. These orderings are constrained by the arguments \code{window} and \code{repeats}, which determine how large the perturbations are and how many are considered. Large values of \code{window} are very time consuming; recommended values are 5 or less, due to the number of permutations which must be considered. Large values of \code{repeats} will eventually converge to an ordering in which all local rearrangements of size \code{window} have been optimized with respect to the number of crossovers. 
#' @return The original object with a new map component. Any pre-existing map will be retained as component $oldmap. 
#' @seealso \code{\link[mpMap]{mpestrf}}, \code{\link[mpMap]{mpgroup}}, \code{\link[mpMap]{add3pt}}, \code{\link[seriation]{seriate}}, \code{\link[qtl]{ripple}}


mporder <-
function(object, chr, type=c("2", "m"), mapfx=c("haldane", "kosambi"), rm.rf=TRUE, window=3, repeats=1, criterion=c("Path_length", "AR_events", "AR_deviations", "Gradient_raw", "Inertia", "Least_squares", "minXO", "lkhdsum"), missfx=2, ...)
{
  if (!inherits(object, "mpcross")) stop("Object must be of class mpcross")

  require(seriation)
  if (missing(object)) 
	stop("Object is required for analysis")

  if (is.null(object$rf))
	stop("Must calculate recombination fractions prior to ordering")

  if (missing(criterion)) criterion <- "Path_length"

  decreasing <- FALSE
  if (criterion %in% c("Path_length", "AR_events", "AR_deviation", "Least_squares", "minXO")) decreasing <- TRUE

  ### what to do if missing values in the recombination fraction matrix
  if (sum(is.na(object$rf$theta))>0 & rm.rf==TRUE) {
	keepmrk <- cleanrf(object)
	obj <- subset(object, markers=keepmrk)
	cat("These markers have been removed due to missing theta estimates: \n")
	cat(setdiff(colnames(object$finals), colnames(obj$finals)), "\n")
	cat("Suggestion is to use add3pt() to insert them into framework map\n")
  } 
  if (rm.rf==FALSE | sum(is.na(object$rf$theta))==0) obj <- object

  if (missing(mapfx)) 	mapfx <- "haldane" 

  if (mapfx=="haldane") mf <- haldaneR2X else mf <- kosambiR2X

  output <- obj
  if (is.null(obj$map) & is.null(obj$lg))
	stop("No grouping of markers input")

  if (is.null(obj$map) & !is.null(obj$lg))
  {
    # create rough obj$map
    obj$map <- list(length=obj$lg$n.groups)
    for (i in 1:obj$lg$n.groups)
    {
	obj$map[[i]] <- rep(0, sum(obj$lg$groups==i, na.rm=TRUE))
	names(obj$map[[i]]) <- colnames(obj$finals)[which(obj$lg$groups==i)]
    }
  }

  if (missing(chr))	chr <- c(1:length(obj$map))

  if (is.character(chr)) chr <- match(chr, names(obj$map))
  # do 2-pt ordering
  if (type=="2")
  {
     order <- list()
     if (criterion=="minXO") {
      write2cross(obj, "tmp", chr=chr)
      cr <- qtl:::readMWril("", "tmp.ril.csv", "tmp.founder.csv", type=attr(obj, "type"))
     }
     newmap <- obj$map 
     for (i in chr)
     {
	 nam <- match(names(obj$map[[i]]), colnames(obj$rf$theta))
	 if (criterion=="lkhdsum") {
		mat1 <- obj$rf$lod[nam, nam] 
		if (sum(is.na(mat1))>0) mat1 <- fill(fill(mat1, missfx), 1)
		diag(mat1) <- 0
		dmat <- as.dist(mat1)
	 }
	 mat <- obj$rf$theta[nam,nam]
	 mat[mat>=.5] <- 0.49
	 if (sum(is.na(mat))>0) mat <- fill(fill(mat, missfx), 1)
	 diag(mat) <- 0
	 if (criterion != "lkhdsum")  dmat <- as.dist(mf(mat))
	
	if (length(obj$map[[i]])>2)
	{
	 ## test out each of the different ordering techniques, pick the 
	 ## one with the shortest path length
	 methods <- c("TSP", "OLO", "ARSA", "Chen", "MDS", "GW", "HC")
	  ser <- sapply(methods, function(x) return(seriate(dmat, method=x))) 
	  o2 <- do.call("rbind", lapply(ser, get_order))
	  o2 <- rbind(1:nrow(mat), o2)
	 if (criterion!="minXO") {
	  criterion1 <- criterion
	  if (criterion1=="lkhdsum") criterion <- "Path_length"
	  crit <- lapply(ser,function(x) return(criterion(dmat,x,criterion)))
	  crit <- c(criterion(dmat, method=criterion), crit)
	  minx <- which.min(unlist(crit))
	  if (!decreasing) minx <- which.max(unlist(crit))
	 } else {
	 crit <- compare_orders(cr, chr=names(obj$map)[[i]], orders=o2, method="countxo")
	 minx <- which.min(crit[,ncol(o2)+1]) }

	 order[[i]] <- o2[minx,1:ncol(o2)]
	} else order[[i]] <- c(1,2)

	mat2 <-  mat[order[[i]], order[[i]]]
	newmap[[i]] <- cumsum(mf(c(0, mat2[row(mat2)==(col(mat2)+1)])))
	names(newmap[[i]]) <- names(obj$map[[i]])[order[[i]]]
     }
     class(newmap) <- "map"
  }

  # do multi-point ordering
  if (type=="m")
  {
    write2cross(obj, "tmp", chr=chr)
    cr <- qtl:::readMWril("", "tmp.ril.csv", "tmp.founder.csv", type=attr(obj, "type"))
    newmap <- list()
    order <- list()
    chr <- match(names(obj$map)[chr], names(cr$geno))
    for (i in chr)
    {
	rip <- ripple(cr, window=window, chr=names(cr$geno)[i])
	nmrk <- nmar(cr)[i]
	cat("Minimum XO for starting order: ", rip[1,nmrk+1], " for best order: ", rip[2,nmrk+1], "\n")
	cr2 <- cr
	order[[i]] <- rip[2, 1:nmrk]
	repeats2 <- repeats

    	while(repeats2 >0 & (rip[1,nmrk+1]!=rip[2,nmrk+1]))
    	{
	  ## actually have to go in and reorder the markers each time
	  cr2$geno[[i]]$data <- cr2$geno[[i]]$data[, order[[i]]]
	  rip <- ripple(cr2, window=window, chr=names(cr$geno)[i])	
	  cat("Minimum XO for starting order: ", rip[1,nmrk+1], " for best order: ", rip[2,nmrk+1],"\n")
	  order[[i]] <- rip[2, 1:nmrk]
	  repeats2 <- repeats2-1
    	}
	## construct new map from ordering
	nam <- match(colnames(cr2$geno[[i]]$data), colnames(obj$rf$theta))
	mat <- obj$rf$theta[nam,nam]
	## what value goes here is important
	mat <- fill(fill(mat, missfx), 1)
	mat[mat==.5] <- .49
	newmap[[i]] <- cumsum(mf(c(0, mat[row(mat)==(col(mat)+1)][1:(length(nam)-1)])))
	names(newmap[[i]]) <- colnames(cr2$geno[[i]]$data)
    }
    names(newmap) <- names(cr$geno)[chr]
    class(newmap) <- "map"
  }
  output$oldmap <- obj$map
  output$map <- newmap
  output <- maporder(output)

  return(output)
}

