library(switchr)
sinfotxt = readLines("./windowsbugdat.txt")


sesstxt <- "R version 3.0.2 (2013-09-25)
Platform: x86_64-w64-mingw32/x64 (64-bit)

locale:
LC_COLLATE=English_United Kingdom.1252
LC_CTYPE=English_United Kingdom.1252
LC_MONETARY=English_United Kingdom.1252
LC_NUMERIC=C
LC_TIME=English_United Kingdom.1252

attached base packages:
stats graphics grDevices utils datasets methods base

other attached packages:
ggplot2_0.9.3.1 R2jags_0.03-12 rjags_3-12 coda_0.16-1
lattice_0.20-29 lme4_1.1-5 Rcpp_0.11.1 Matrix_1.1-2
reshape2_1.2.2 knitr_1.5

loaded via a namespace (and not attached):
abind_1.4-0 boot_1.3-9 colorspace_1.2-4
dichromat_2.0-0 digest_0.6.4 evaluate_0.5.1
formatR_0.10 grid_3.0.2 gtable_0.1.2
highr_0.3 labeling_0.2 MASS_7.3-29
minqa_1.2.3 munsell_0.4.2 nlme_3.1-111
parallel_3.0.2 plyr_1.8.1 proto_0.3-10
R2WinBUGS_2.1-19 RColorBrewer_1.0-5 RcppEigen_0.3.2.0.2
scales_0.2.3 splines_3.0.2 stringr_0.6.2
tools_3.0.2"

# This needs to be formatted to match txt above
sesstxt <- strsplit(sesstxt, split = '\n')[[1]]

removeLib("winbugtst")
switchTo("winbugtst", seed = sesstxt)



switchTo("emptytest")
man = makeSeedMan(parseSessionInfoString(sesstxt))
repourl = lazyRepo(man)

availpkgs = available.packages(repourl, type="source")
dim(availpkgs)
pkgs = availpkgs[,"Package"]
pkgs

install_packages(pkgs, repourl)

installed.packages()
