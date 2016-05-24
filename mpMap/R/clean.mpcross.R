#' Check data format and compute summary statistics for genotypes
#' 
#' Given an object of class 'mpcross', the function checks that the data is in the correct format, containing founder and final genotypes, ids, and a pedigree. The number of markers genotyped for both founders and finals should coincide. The pedigree should be completely numeric. Markers which are not polymorphic across the founders are removed, as are markers which have missing values in the founders. 
#'
#' Summary statistics for the genotypes are printed, included the number of markers with varying levels of missing data, with varying levels of segregation distortion, and with different numbers of alleles. 
#' @export 
#' @import qtl
#' @S3method clean mpcross
#' @method clean mpcross
#' @aliases mpsegrat
#' @param object Object of class \code{mpcross}
#' @param ... Additional arguments
#' @return 
#' \item{alleles}{Number of alleles at each marker}
#' \item{missing}{Percent missing data at each marker}
#' \item{seg}{Matrix with one row for each marker and columns for the marker name, the chisquare test for segregation distortion, and the p-value of the test}
#' @seealso \code{\link[mpMap]{mpcross}}
#' @examples
#' map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, qtl=matrix(data=c(1, 45, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE),seed=1)
#' dat.chk <- clean(sim.dat)

clean.mpcross <- function(object, ...)
{

 if (is.null(object$founders))
	stop("Founder genotypes missing in object") 

 if (is.null(object$finals))
	stop("Final genotypes missing in object")

 if (is.null(object$pedigree))
	stop("Pedigree missing in object")

 if (is.null(object$id))
	stop("Genotype IDs missing in object") 

 n.founders <- nrow(object$founders)
 n.mrk <- ncol(object$founders)
 
 if (n.mrk != ncol(object$finals))
	stop("Number of markers for finals and founders does not match")

 # Check pedigree format
 pn <- check_ped(object$pedigree)

 # Remove markers which do not differ among founders
 fdr.alleles <- apply(object$founders, 2, function(x) return(length(table(x))))
 drop1 <- which(fdr.alleles==1)

 # Remove markers which have missing founder values
 fdr.missing <- apply(object$founders, 2, function(x) return(sum(is.na(x))))
 drop2 <- which(fdr.missing>0)

 drop <- unique(c(drop1, drop2))
 if (length(drop)>0)
 object <- subset(object, setdiff(colnames(object$finals), drop))

 # Print off summary statistics
	# no. biallelic markers, no. multiallelic markers
	# no. markers with >5%, >10%, >20% missing data
 pctmiss <- apply(object$finals, 2, function(x) return(sum(is.na(x))/length(x)))
	# no. markers with segregation p-values <1e-5, <1e-10, etc.
 segrat <- mpsegrat(object)
 segpval <- segrat[,3]

 return(list(drop1=drop1, drop2=drop2, alleles=fdr.alleles, missing=pctmiss, seg=segrat))
}


