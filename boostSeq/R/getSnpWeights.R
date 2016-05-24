getSnpWeights <- function(snps, gwas.resultfile, pColname = "P", snpColname = "SNP") {
  p.dat  <- read.table(gwas.resultfile, header = TRUE, stringsAsFactors = FALSE)
  p.sel  <- p.dat[p.dat[, snpColname] %in% snps, c(snpColname, pColname)]
  p.log  <- log10(p.sel[, pColname])
  names(p.log) <- p.sel[, snpColname]
  return(abs(p.log - min(p.log) +1))
}


# This function extracts the risk-associated allele from PLINK result files (.tdt and .assoc files)
# returns a vector of preferred nucleotides for each SNP, named by SNP IDs
getRefAlleles <- function(file) {
  cat("Extracting reference alleles...\n")
  alleles.pref <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  alleles.ref        <- c(
    alleles.pref[alleles.pref$OR > 1, "A1"], 
    alleles.pref[alleles.pref$OR <= 1, "A2"]
    )
  names(alleles.ref) <- c(
    alleles.pref[alleles.pref$OR > 1, "SNP"],
    alleles.pref[alleles.pref$OR <= 1, "SNP"]
    )
  return(alleles.ref)
}
