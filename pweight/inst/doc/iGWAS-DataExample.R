## ------------------------------------------------------------------------
#uncomment on Windows
#setInternet2(use = TRUE)
download.file("http://github.com/dobriban/pweight/raw/master/vignettes/iGWAS.cad.example.RData", destfile = "exampledata", method = "libcurl")
load("exampledata")

## ------------------------------------------------------------------------
c4d_ss <- 4/(1/c4d_ncases + 1/c4d_ncontrols)
cardiogram_ss <- 4/(1/cardiogram_ncases + 1/cardiogram_ncontrols)

## ------------------------------------------------------------------------
source("../R/iGWAS.R")
source("../R/bayes_weights.R")
res_default <- iGWAS(P_current=cad$P_c4d, N_current=c4d_ss, P_prior=cad$P_cardiogram, N_prior=cardiogram_ss)

## ------------------------------------------------------------------------
cat(c("Number of genome-wide significant SNPs (P < 5e-8) with no weighting: ", sum(cad$P_c4d < 5e-8) ))

## ------------------------------------------------------------------------
w = res_default$w
plot(-log10(cad$P_cardiogram), w, main="How Bayes weights depend on p-values", ylab=expression(italic(w)), xlab=expression(-log[10](italic(p))))
cat(c("Maximum weight reached at P-value of ", cad$P_cardiogram[which.max(w)] ))

## ------------------------------------------------------------------------
res_default <- iGWAS(P_current=cad$P_c4d, N_current=c4d_ss, P_prior=cad$P_cardiogram, N_prior=cardiogram_ss,figure="T", GWAS_data_frame=cad)

## ------------------------------------------------------------------------
source("../R/spjotvoll_weights.R")
phi=seq(0,1,by=0.1)
numSNP <- sapply(phi,function(x){sum(iGWAS(P_current=cad$P_c4d, N_current=c4d_ss, P_prior=cad$P_cardiogram, N_prior=cardiogram_ss, phi=x)$sig_ind)})
plot(phi, numSNP, type="o", main="Choice of phi impacts the number of significant SNPs", ylab="# SNPs", xlab=expression(phi))
cat(c("Choice of phi that maximizes the number of genome-wide significant SNPs:", phi[which.max(numSNP)] ))

