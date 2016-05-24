
set.genomic.sex <- function(x, chr.x = 23, chr.y = 24, plot = FALSE, verbose = getOption("gaston.verbose",TRUE)) {
  if( !all(c("hz", "maf", "callrate") %in% names(x@snps) )) {
    if(verbose) cat("Computing basic stats\n", set.ped_stats = FALSE, set.p = FALSE, set.mu_sigma = FALSE)
    x <- set.stats(x)
  }
  if(verbose) cat("Extracting SNPs from chr X and chr Y SNPs\n")
  X <- x[, x@snps$chr == chr.x & !is.na(x@snps$maf) & x@snps$maf > 0]
  Y <- x[, x@snps$chr == chr.y & !is.na(x@snps$maf) & x@snps$maf > 0]
  if(!("callrate" %in% names(X@ped)))
    X <- set.stats(X, verbose = FALSE, set.snps_stats = FALSE, set.p = FALSE, set.mu_sigma = FALSE)
  if(!("hz" %in% names(Y@ped)))
    Y <- set.stats(Y, verbose = FALSE, set.snps_stats = FALSE, set.p = FALSE, set.mu_sigma = FALSE)

  y.callrate <- Y@ped$callrate
  x.hz <- X@ped$hz
  x@ped$genomic.sex <- kmeans( x = cbind(y.callrate, x.hz), centers = matrix( c(max(y.callrate), min(x.hz), min(y.callrate), max(x.hz)), ncol = 2 ) )$cluster
  if(plot) {
    plot( x.hz, y.callrate, xlab = "X heterozygosity", ylab = "Y callrate", col = x@ped$genomic.sex, main = "Determination of Genomic Sex" )
    legend("bottomleft", pch = 1, col = 1:2, c("M","F"))
  }
  if(verbose) {
    cat("Slot @ped$genomic.sex has been set\n")
    cat("There are", sum(!is.na(x@ped$sex) & x@ped$sex == 1 & x@ped$genomic.sex == 2), "individuals with sex = M and genomic sex = F\n")
    cat("There are", sum(!is.na(x@ped$sex) & x@ped$sex == 2 & x@ped$genomic.sex == 1), "individuals with sex = F and genomic sex = M\n")
  }
  x
}
