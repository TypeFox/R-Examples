#' Construct mpcross objects from datafiles
#'
#' Generate an mpcross object by reading in components from files - requires founders marker data, finals marker data, pedigree, and IDs for all lines. Marker map and phenotypic data are optional. 
#' @export
#' @param founderfile File containing founder genotypes - should have 1+(number of founders) columns. The first column contains the marker names - first space left blank. The first row for the other columns contains the founder name. Additional rows contain observed marker data for all founders.
#' @param finalfile File containing final genotypes - should have 1+(number of lines) columns. The first column contains the marker names - first space left blank. The first row for the other columns contains line names. Additional rows contain observed marker data for all lines.
#' @param pedfile File containing pedigree data - should have three or four columns - first three columns indicate id, mother and father; fourth column is a flag for whether the lines was observed. 
#' @param mapfile File containing linkage map - should contain three columns - one for marker names, one for chromosome assignment and one for map position in cM
#' @param phenofile File containing phenotypic data - should contain one column for each phenotype, with rows indexing lines.
#' @seealso \code{\link[mpMap]{sim.mpcross}}, \code{\link[mpMap]{mpcross}}

read.mpcross <- function(founderfile, finalfile, pedfile, mapfile, phenofile)
{
  ## what format? 
  founders <- read.table(founderfile)
  finals <- read.table(finalfile)

	if(!missing(pedfile))
	{
		ped1 <- read.table(pedfile)
		## potentially need to check.ped as well. 
		ped <- convertped(ped1)
	}
	else ped <- NULL

  id <- match(rownames(finals), ped1[,1])
  if (ncol(ped)==4) id <- which(ped[,4]==1)
#  fid <- rownames(founders) 
  fid <- 1:nrow(founders)

  object <- mpcross(founders, finals, ped, id, fid)
  
  if (!missing(mapfile))
  {
   	mapin <- read.table(mapfile, header=TRUE)
  	if (ncol(mapin)!=3) stop("Map file is incorrectly formatted")

  	map <- list()
  	for (i in names(table(mapin[,2])))
  	{
	  map[[i]] <- mapin[which(mapin[,2]==i), 3]
	  names(map[[i]]) <- mapin[which(mapin[,2]==i), 1]
  	}
  	object$map <- map
  }

  if (!missing(phenofile))
  {
	phein <- read.table(phenofile, header=TRUE)
  	m <- match(rownames(object$finals), rownames(phein))
	if (sum(is.na(m))>0) cat("Lines in finals data which are not in phenotypes")
	if (length(m) < nrow(phein)) cat("Lines in phenotypes which are not in finals data")
	m <- m[!is.na(m)]
	
	pheno <- phein[m,]
	object$pheno <- as.data.frame(pheno)
  }

  object
}

