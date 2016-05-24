## Functions dealing with setting up and ordering genotypes

enumerate_genotypes <- function(nloci = 1, locinames = NULL, allelenames = NULL) {

  ## insert default loci and allele names if NULL, expand allelenames 
  if (is.null(locinames)) {
    locinames <- paste('Locus_', 1:nloci, sep = "") 
  }
  if (is.null(allelenames)) { 
    allelenames <- rep(c("1", "2"), nloci)
    dim(allelenames) <- c(2, nloci)
    allelenames <- t(allelenames)
  }
  if (length(allelenames) == 2 & nloci > 1) {
    allelenames <- rep(allelenames, nloci)
    dim(allelenames) <- c(2, nloci)
    allelenames <- t(allelenames)
  }

  # check consistency
  if (nloci != length(locinames)) { 
    stop('length(locinames) must be equal to nloci.') 
  }
  if (nloci != dim(allelenames)[1]) { 
    stop('dim(allelenames)[1] must be equal to nloci.') 
  }

  # construct data frame with 3^nloci genotypes
  tmpargs <- NULL 
  for (i in 1:nloci) {
    tmpargs[[paste('Locus', i, sep = "")]] <- 
      factor(c(paste(allelenames[i, c(1, 1)], collapse = ""), 
               paste(allelenames[i, c(1, 2)], collapse = ""), 
               paste(allelenames[i, c(2, 2)], collapse = "")))
  } 
  genotypes <- do.call(expand.grid, list(tmpargs,KEEP.OUT.ATTRS=FALSE), quote=TRUE)
  colnames(genotypes) <- locinames
  return(genotypes)
}

partial_genotype_order <- function(plusallele) {
  ## build up partial ordering of genotype space based on input vector plusallele

  
  apart <- NULL  
  nloci <- length(plusallele)            
    for (locus in 1:nloci) {
    g <- 1:(3^nloci) 
    dim(g) <- rep(3, nloci) 
    sortg <- aaply(g, locus, sort) 
    dim(sortg) <- c(3 ^ (nloci - 1), 3)
    if (plusallele[locus] == 2) {
      apart <- rbind(apart, cbind(sortg[, 1], sortg[, 2]), 
                     cbind(sortg[, 2], sortg[, 3]))
    }
    if (plusallele[locus] == 1) {
      apart <- rbind(apart, cbind(sortg[, 3], sortg[, 2]), 
                     cbind(sortg[, 2], sortg[, 1]))
    }
  }
  return(apart)
}
