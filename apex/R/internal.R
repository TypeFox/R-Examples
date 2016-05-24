
##
## Internal functions: not documented, not exported
##
## Thibaut Jombart, April 2015
##

## compute the number of missing sequences for each gene in a list of DNAbin matrices or phyDat objects 'x'
.nMissingSequences <- function(x){
    ## only keep non-empty matrices
    x <- lapply(x, as.character)
    x <- x[sapply(x, nrow)>0]
    out <- sapply(x, function(e) sum(apply(e=="-",1,all)))
    return(sum(out,na.rm=TRUE))
}


# Compare requested locus names with actual names in multidna object
.checkLocusNames <- function(x, loci = NULL) {
  # check that locus names can be found
  if(!is.null(loci)) {
    missing <- setdiff(loci, getLocusNames(x))
    if(length(missing) > 0) {
      missing <- paste(missing, collapse = ", ")
      warning(paste("The following loci could not be found:", missing), call. = FALSE)
    }
    loci <- intersect(loci, getLocusNames(x))
  } else if(is.logical(loci) | is.numeric(loci)) {
    loci <- getLocusNames(x)[loci]
  }
  # set to all locus names if none specified
  if(is.null(loci)) loci <- getLocusNames(x)
  # return NULL if no locus names match
  if(length(loci) == 0) return(NULL)
  loci
}

# Compare requested ids with actual ids in DNAbin object
.checkIDs <- function(dna, ids = NULL) {
  # check that id names can be found
  if(!is.null(ids)) {
    missing <- setdiff(ids, labels(dna))
    if(length(missing) > 0) {
      missing <- paste(missing, collapse = ", ")
      warning(paste("The following ids could not be found:", missing), call. = FALSE)
    }
    ids <- intersect(ids, labels(dna))
  } else if(is.logical(ids) | is.numeric(ids)) {
    ids <- labels(dna)[ids]
  }
  # set to all locus names if none specified
  if(is.null(ids)) ids <- labels(dna)
  # return NULL if no locus names match
  if(length(ids) == 0) return(NULL)
  ids
}


# Return a vector of logicals denoting if each DNAbin sequence contains only gaps
.isGapOnly <- function(dna) {
  sapply(as.character(as.list(dna)), function(this.seq) all(this.seq == "-"))
}
# the same for phyDat objects
.isGapOnlyPhyDat <- function(dna, gap="-") {
    allLevels <- attr(dna, "allLevels")
    ind <- match(gap, allLevels)
    sapply(dna, function(this.seq) all(this.seq == ind))
}
