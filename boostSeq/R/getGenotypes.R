getGenotypes <- function(
  snps, 
  ped = NULL, 
  map = NULL, 
  remove.homozygous = FALSE, 
  as.geno.objects = TRUE) 
{
  
  if(missing(snps))
    stop("Parameter 'snps' has to be specified")
  
  snps <- as.vector(snps)
  
  if(is.null(snps) || length(snps) < 1) {
    warning("No SNPs specified to retrieve genotypes for")
    return(NULL)
  }
  
  if(is.null(ped) || is.null(map)) {
    stop("Please specify ped/map files for genotype retrieval.\n")
  }
  
  gts <- NULL
  
  cat("Reading genotype files ...\n")
  mapinfo <- readMapfile(map)
  snp.colmask <- mapinfo$SNP %in% snps
  
  # read ped line by line, extract SNPs of interest
  pedfile <- file(ped, open = "r")
  ped.line <- readLines(pedfile, n = 1)
  
  while( length(ped.line) != 0 ) {
    ped.line.vec <- unlist(strsplit(ped.line, split = "\\s"))
    # remove leading non-genotype columns
    gts.line <- ped.line.vec[(length(ped.line.vec) - 2*length(snp.colmask) +1):length(ped.line.vec)]
    # select target snps and combine the two alleles for each SNP
    gts.line.selected <- sapply(which(snp.colmask) *2, function(idx) paste(gts.line[idx -1], gts.line[idx]))
    # format for 'genotype' class (recode numeric nucleotides to character)
    gts.line.selected <- sub(" ", "", gts.line.selected, fixed = TRUE)
    gts.line.selected <- gsub("1", "A", gts.line.selected, fixed = TRUE)
    gts.line.selected <- gsub("2", "C", gts.line.selected, fixed = TRUE)
    gts.line.selected <- gsub("3", "G", gts.line.selected, fixed = TRUE)
    gts.line.selected <- gsub("4", "T", gts.line.selected, fixed = TRUE)
    # add line to whole gts matrix and read next line from file
    gts <- rbind(gts, gts.line.selected)
    ped.line <- readLines(pedfile, n = 1)
  }
  close(pedfile)
  colnames(gts) <- mapinfo$SNP[mapinfo$SNP %in% snps]
  
  # replace everything except nucleotides with NA
  gts[grep("[^ACGT]", gts)] <- NA
  
  if( remove.homozygous ) {
    # tabulate all gts gives number of homo- and heterozygous - so we just extract those that have more than one 'table' entry
    hom.het.hom.count <- apply(gts, 2, table)
    # apply returns matrix if all are same zygosity state (same number of 'table' entries), convert to list in this case
    if( class(hom.het.hom.count) == "matrix" ) 
      hom.het.hom.count <- as.data.frame(hom.het.hom.count)
    gts.zygosity <- lapply(hom.het.hom.count, length)
    gts <- gts[, unlist(gts.zygosity) != 1]
  }
  
  if( as.geno.objects ) {
    if(is.null(gts))
      return(NULL)
    if( remove.homozygous ) 
      if( class(hom.het.hom.count) != "matrix" )
        gts <- as.matrix(gts)
    # remove all SNPs that have only NA genotypes (converting to class genotype fails in this case)
    gts.cols.nacount <- apply(is.na(gts),2,sum)
    gts <- gts[, !(gts.cols.nacount == nrow(gts))]
    cat("Creating genotype objects...\n")
    geno <- lapply(as.data.frame(gts), function(x) { genotype(x, sep = "") })
    if(length(geno) == 1) {
      # genotype conversion by data frame does not retain names for a single element
      names(geno) <- mapinfo$SNP[mapinfo$SNP %in% snps]
    }
    return(geno)
  } else {
    return(gts)
  }
}

readMapfile <- function(mapfile) {
  mapinfo <- read.table(mapfile, header = FALSE, row.names = NULL)
  if(ncol(mapinfo) == 4) mapinfo <- mapinfo[, c(1,2,4)]
  if(ncol(mapinfo) != 3) 
    stop("Mapfile format unknown")
  colnames(mapinfo) <- c("CHR", "SNP", "BP")
  return(mapinfo)
}

