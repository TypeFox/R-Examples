library("ACNE")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DATA: Lx2xI allele-specific signals for six different SNPs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
filenames <- sprintf("V%d.Rbin", 1:6)
pathnames <- system.file("extData", filenames, package="ACNE")
Ys <- lapply(pathnames, FUN=function(p) snpMatrixToArray(loadToEnv(p)$V))
names(Ys) <- sprintf("SNP #%d", seq_along(Ys))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ACNE fitting of NMF to the six SNPs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (name in names(Ys)) {
  Y <- Ys[[name]]
  fit <- fitSnpNmfArray(Y)
  str(fit)
  plot(fit, lim=c(0,2^14), main=name)
}
