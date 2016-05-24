#' Add markers onto an existing 'mpcross' object
#' 
#' If the 'mpcross' object contains only genotypes and pedigree, the markers are added to the genotypes after some format checking. If the 'mpcross' object contains estimated recombination fractions and LOD scores, the markers are added in after estimating recombination fractions with the existing markers. If the 'mpcross' object contains a map, the markers are added onto the map using 3-point mapping. 
#' @export
#' @param mpcross1 Original object of class \code{mpcross}
#' @param mpcross2 Additional object of class \code{mpcross}
#' @param r Grid of recombination fractions over which to maximize likelihood
#' @param theta Threshold for recombination fractions used in constructing linkage groups
#' @param LOD Threshold for LOD scores used in constructing linkage groups
#' @param mapfx Map function for converting recombination distances to map distances
#' @return Object of class 'mpcross' which is of the same stage as the first object. Recombination fractions and linkage groups will be recomputed after adding in the markers in the second object if necessary. 
#' @seealso \code{\link[mpMap]{mpestrf}}, \code{\link[mpMap]{mpgroup}}, \code{\link[mpMap]{mpcross}}

mpadd <-
function(mpcross1, mpcross2, r, theta=.15, LOD=5, mapfx=c("haldane", "kosambi"))
{
  output <- list()

  if (missing(mapfx)) mapfx <- "haldane"

  if (mapfx=="kosambi") mf <- kosambiR2X else mf <- haldaneR2X

  if (missing(mpcross1) | missing(mpcross2))
	stop("Missing a required argument for this function")

  if (missing(r)) r <- c(0:20/200, 11:50/100)

  # Check for overlapping markers
  mrk.isect <- intersect(colnames(mpcross1$founders), colnames(mpcross2$founders))
  data2 <- subset(mpcross2, setdiff(colnames(mpcross2$finals), mrk.isect))

  # Check that the individuals are the same and in the same order
  if (any(sort(mpcross1$id)!=sort(data2$id)))  
	stop("Different observations in the two objects")
  
  if(any(mpcross1$id!=data2$id)) {
	data2$finals <- data2$finals[match(data2$id, mpcross1$id),]
 	data2$id <- mpcross1$id
	}

  # Check that the founders are the same and in the same order
  if (any(sort(mpcross1$fid)!=sort(data2$fid)))
	stop("Different founders in the two objects")
 
  if (any(mpcross1$fid!=data2$fid)) {
	data2$founders <- data2$founders[match(data2$fid, mpcross1$fid),]
	data2$fid <- mpcross1$fid 
	}

  # Check that the two objects have the same pedigree
  if (any(mpcross1$ped != data2$ped))
	stop("Different pedigrees in the two objects")

  # start by combining the founders and finals
  output$founders <- cbind(mpcross1$founders, data2$founders)
  output$finals <- cbind(mpcross1$finals, data2$finals)
  output$pedigree <- mpcross1$pedigree
  output$id <- mpcross1$id
  output$fid <- mpcross1$fid

  nmrk1 <- ncol(mpcross1$finals)
  nmrk2 <- ncol(mpcross2$finals)

  # combine recombination fraction matrices
  if (!is.null(mpcross1$rf))  
  {
    if(is.null(mpcross2$rf))
	data2 <- mpestrf(data2, r) 
    output$rf <- combine_rf(mpcross1, data2, r)
    if (!is.null(mpcross1$ld))
    {
      if(is.null(mpcross2$ld))
	data2 <- mpcalcld(data2)
      output$ld <- combine_ld(mpcross1, data2, output$rf$theta[1:nmrk1, nmrk1+1:nmrk2])
    }
  }

  if (!is.null(mpcross1$lg))	output$lg1 <- mpcross1$lg
  if (!is.null(mpcross2$lg))	output$lg2 <- mpcross2$lg
  if (!is.null(mpcross1$map))   output$map1 <- mpcross1$map
  if (!is.null(mpcross2$map)) 	output$map2 <- mpcross2$map 

  class(output) <- "mpcross"
  return(output)
}

