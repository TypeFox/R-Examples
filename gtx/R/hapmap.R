# a default of path = "ftp://ftp.ncbi.nlm.nih.gov/hapmap"
# would work if read.table(gzfile("ftp://...")) worked

hapmap.read.haplotypes <- function (chr, path, sample = "CEU", release = 22) {
  if (release == 21) {

    legend <- read.table(gzfile(paste(path, "/phasing/2006-07_phaseII/phased/genotypes_chr", chr, "_", sample, "_r21_nr_fwd_legend.txt.gz", sep = "")), comment.char = "", as.is = TRUE, header = TRUE)
    return(list(legend = legend, haplotypes = matrix(scan(gzfile(paste(path, "/phasing/2006-07_phaseII/phased/genotypes_chr", chr, "_", sample, "_r21_nr_fwd_phased.gz", sep = "")), what = integer(0), quiet = TRUE), nrow = nrow(legend))))
    
  } else if (release == 22) {

    legend <- read.table(gzfile(paste(path, "/phasing/2007-08_rel22/phased/genotypes_chr", chr, "_", sample, "_r22_nr.b36_fwd_legend.txt.gz", sep = "")), comment.char = "", as.is = TRUE, header = TRUE)
    return(list(legend = legend, haplotypes = matrix(scan(gzfile(paste(path, "/phasing/2007-08_rel22/phased/genotypes_chr", chr, "_", sample, "_r22_nr.b36_fwd.phase.gz", sep = "")), what = integer(0), quiet = TRUE), nrow = nrow(legend))))

  } else {
    stop("release must be 21 or 22")
  }
}

## hapmap22 <- lapply(1:22, hapmap.read.haplotypes, path = "~/data/hapmap")

hapmap.snpdata <- function(params, hapmap) {
  snpinfo <- data.frame(NULL)
  nindiv <- sapply(1:22, function(chr) return(ncol(hapmap[[chr]]$haplotypes)))
  stopifnot(all(nindiv == nindiv[1]))
  snpdata <- list(snpinfo = data.frame(NULL), data = data.frame(HID = 1:nindiv[1]), ploidy = 1)  
  for (chr in 1:22) {
    hm <- match(params$snp, hapmap[[chr]]$legend$rs)
    hmi <- na.omit(hm)
    if (length(hmi) > 0) {
      snpdata$snpinfo <- rbind(snpdata$snpinfo, hapmap[[chr]]$legend[hmi, c("rs", "X1", "X0")])
      hdata <- as.data.frame(t(hapmap[[chr]]$haplotypes[hmi, , drop = FALSE]))
      names(hdata) <- paste(hapmap[[chr]]$legend$rs[hmi], hapmap[[chr]]$legend$X1[hmi], sep = "_")
      snpdata$data <- cbind(snpdata$data, hdata)
    }
  }
  names(snpdata$snpinfo) <- c("snp", "coded.allele", "noncoded.allele")
  stopifnot(all(paste(snpdata$snpinfo$snp, snpdata$snpinfo$coded.allele, sep = "_") == names(snpdata$data)[-1]))
  snpdata$snpinfo$coded.freq <- apply(snpdata$data[ , -1], 2, mean)

  params$codeB <- paste(params$snp, params$coded.allele, sep = "_")
  if(!all(params$codeB %in% names(snpdata$data))) snpdata <- align.snpdata.coding(params, snpdata, ploidy = 1)$snpdata
  stopifnot(all(params$codeB %in% names(snpdata$data))) # this should not happen
  snpdata$data <- subset(snpdata$data, select = params$codeB)
  return(snpdata)
}
