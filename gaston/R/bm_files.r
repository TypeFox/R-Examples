read.bed.matrix <- function(basename, bed = paste(basename, ".bed", sep=""), fam = paste(basename, ".fam", sep=""), 
                                      bim = paste(basename, ".bim", sep=""), rds = paste(basename, ".rds", sep=""),
                            verbose = getOption("gaston.verbose",TRUE)) {

  bed <- path.expand(bed)
  if(!file.exists(bed)) { # peut-être on a donné le .bed pour basename
    if(grep("\\.bed$", basename)) {
      basename <- sub("\\.bed$", "", basename)
      bim <- paste(basename, ".bim", sep="")
      fam <- paste(basename, ".fam", sep="")
      bed <- paste(basename, ".bed", sep="") 
      if(!is.null(rds)) rds <- paste(basename, ".rds", sep="")
    }
  }

  if(!file.exists(bed)) stop("file ", bed, " not found")

  if(is.null(rds) || !file.exists(rds)) {
    if(!file.exists(fam)) stop("file ", fam, " not found")
    if(!file.exists(bim)) stop("file ", bim, " not found")

    if(verbose) cat("Reading", fam, "\n")
    ped <- read.table(fam, stringsAsFactors = FALSE) 
    colnames(ped) <- pednames

    if(verbose) cat("Reading", bim, "\n")
    snp <- read.table(bim, stringsAsFactors = FALSE)
    colnames(snp) <- snpnames

    if(verbose) cat("Reading", bed, "\n")
    bed <- .Call('gg_read_bed_file', bed, nrow(ped), nrow(snp))
    x <- new("bed.matrix", bed = bed, snps = snp, ped = ped,                            
      p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
      standardize_mu_sigma = FALSE )
    if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
    return(x)
  }
 
  # reading rds [on ne fait pas set.stats dans ce cas !!]
  if(verbose) cat("Reading", rds, "\n")
  x <- readRDS(rds)
  if ( is(x) != "bed.matrix" ) stop("The object in file ", rds, " is not a bed.matrix")

  if(verbose) cat("Reading", bed, "\n")
  x@bed <- .Call('gg_read_bed_file', bed, nrow(x@ped), nrow(x@snps))

  x
}


write.bed.matrix <- function(x, basename, bed = paste(basename, ".bed", sep=""), fam = paste(basename, ".fam", sep=""),
                                          bim = paste(basename, ".bim", sep=""), rds = paste(basename, ".rds", sep="")) {
  if ( is(x) != "bed.matrix" ) stop("x must be a bed.matrix")

  if(!is.null(fam)) {
    if(all(pednames %in% names(x@ped)))
      write.table(x@ped[,pednames], file = fam, row.names=FALSE, col.names=FALSE, quote = FALSE)
    else
      warning("Can't create .fam file from x: missing columns in x@ped")
  }

  if(!is.null(bim)){
    if(all(snpnames %in% names(x@snps)))
      write.table(x@snps[,snpnames], file = bim, row.names=FALSE, col.names=FALSE, sep = "\t", quote = FALSE)
    else
      warning("Can't create .bim file from x: missing columns in x@snps")
  }

  if(!is.null(rds))
    saveRDS(x, rds)

  bed <- path.expand(bed)
  invisible(.Call('gg_write_bed_file', x@bed, bed))
}




read.vcf <- function(filename, max.snps, verbose = getOption("gaston.verbose",TRUE)) {
  filename <- path.expand(filename)
  xx <- vcf_open(filename)
  if(is.null(xx)) stop("File not found")
  samples <- vcf_getsamples(xx) 
  vcf_selectsamples( xx, samples )

  d <- ""
  f <- function() vcf_readLineRaw(xx, d) 

  if(missing(max.snps)) {
    if(verbose) cat("Counting diallelic variants\n")
    max.snps <- .Call("gg_count_dia_vcf", PACKAGE="gaston", f, d)
    vcf_close(xx);
    vcf_reopen(xx);
  }

  if(verbose) cat("Reading", max.snps, "diallelic variants for", length(samples), "individuals\n");
  
  L <- .Call("gg_read_vcf", PACKAGE="gaston", f, d, length(samples), max.snps)
  vcf_close(xx)
  snp <- data.frame(chr = L$chr, id = L$id, dist = 0, pos = L$pos , A1 = L$A1, A2 = L$A2, stringsAsFactors = FALSE)
  ped <- data.frame(famid = samples, id = samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}

count.vcf <- function(filename) {
  xx <- vcf_open(filename)
  if(is.null(xx)) stop("File not found")
  samples <- vcf_getsamples(xx) 
  vcf_selectsamples( xx, samples )
  d <- ""
  f <- function() vcf_readLineRaw(xx, d) 
  n <- .Call("gg_count_dia_vcf", PACKAGE="gaston", f, d)
  vcf_close(xx)
  n
}


