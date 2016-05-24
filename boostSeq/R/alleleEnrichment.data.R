# pheno data is appended to the bottom of the matrix
alleleEnrichment.data.alleles <- function(
  refalleles = NULL,
  map, 
  ped,
  pheno.file = NULL, 
  ped.sampleid.column,
  na.samples = na.omit,
  na.snps = na.pass
) {
  
  # read genotypes
  # mindest-allelfrequenz? was ist mit nur homozygoten?
  alleles <- getGenotypes(snps = names(refalleles), ped = ped, map = map, as.geno.objects = T)
  
  if(length(alleles) != length(refalleles))
    stop(paste("Could not read all SNPs from pedfile:", paste(names(refalleles)[!(names(refalleles) %in% names(alleles))], collapse = ";")))
  
  if(is.null(refalleles)) {
    # construct zygosity matrix
    # each snp get three rows: (hom/het/hom coded with 0 and 1)
    alleles <- t(as.data.frame(sapply(
                 names(alleles), 
                 function(snp) {
                   cbind(
                     as.numeric(allele.count(alleles[[snp]])[, 1] == 0), 
                     as.numeric(allele.count(alleles[[snp]])[, 1] == 1), 
                     as.numeric(allele.count(alleles[[snp]])[, 1] == 2)
                   )
                 },
                 simplify = FALSE
               )))
  
  } else {
  
    # allele count matrix express alleles matrix in terms of ref allele count, transpose so that samples = columns
    alleles <- t(as.data.frame(sapply(
               names(alleles), 
               function(snp) return(allele.count(alleles[[snp]], refalleles[[snp]])),
               simplify = FALSE
             )))
  }

  # set colnames (sample IDs)
  # read sample ids from pedfile (line by line)
  cat("Extracting sample IDs from pedfile...\n")
  samples <- NULL
  pedfile <- file(ped, open = "r")
  ped.line <- readLines(pedfile, n = 1)
  while( length(ped.line) != 0 ) {
    ped.line.vec <- unlist(strsplit(ped.line, split = "\\s"))
    samples <-c(samples, ped.line.vec[ped.sampleid.column])
    ped.line <- readLines(pedfile, n = 1)
  }
  close(pedfile)
  colnames(alleles) <- samples

  if(!is.null(pheno.file)) {
    phenodat <- alleleEnrichment.data.pheno(pheno.file, ped.sampleid.column)
    # remember rownames because ordering on a one-row matrix coverts to unnamed vector (auto-cast of types in R... hmpf)
    phen.rn <- colnames(phenodat)
    phenodat <- t(as.matrix(phenodat[order(rownames(phenodat)), ]))
    snp.rn <- rownames(alleles) 
    alleles <- alleles[, order(colnames(alleles))]
    alleles <- rbind(alleles, phenodat)
    rownames(alleles) <- c(snp.rn, phen.rn)
  }

  # remove NAs
  # handle rows that contain NAs snp-wise
  alleles <- na.snps(alleles)
  # handle rows that contain NAs sample-wise
  alleles <- t(na.samples(t(alleles)))
  # handle rows that contain NAs genotype-wise
  if(any(is.na(alleles)))
    alleles[is.na(alleles)] <- round(runif(sum(is.na(alleles)))) # paste randomly 0 or 1
  
  return(alleles)
  
}



alleleEnrichment.data.pheno <- function(file, idcol) {
  cat("Extracting phenotype data...\n")
  pheno <- read.table(file, stringsAsFactors = FALSE)
  cn <- pheno[1, (idcol+1):ncol(pheno)]
  pheno <- pheno[-1, ]
  rn <- pheno[, idcol]
  # as.data.frame converts a single column to matrix and not vector
  pheno <- as.data.frame(pheno[, (idcol+1):ncol(pheno)])
  tryCatch( 
    {pheno <- as.data.frame(lapply(pheno, as.numeric))},
    warning = function(w) stop("Covariates contain values that cannot be converted to numeric.\n")
  )
  colnames(pheno) <- cn
  rownames(pheno) <- rn
  return(pheno)
}

