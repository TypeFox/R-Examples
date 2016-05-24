## ---- echo = FALSE-------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", cache = FALSE,
                      fig.width = 5, fig.height = 4) # warning = FALSE

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("inbeedR")

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("mastoffel/inbreedR", build_vignettes = TRUE)

## ------------------------------------------------------------------------
library("inbreedR")

## ---- eval = FALSE-------------------------------------------------------
#  ?inbreedR

## ---- echo=TRUE, results='hide'------------------------------------------
library(inbreedR)
data("mouse_msats") # microsatellite data 
data("mouse_snps")  # snp data
data("bodyweight")  # fitness data


## ------------------------------------------------------------------------
data("mouse_snps")
mouse_snps[1:10, 1:10]

## ------------------------------------------------------------------------
check_data(mouse_snps, num_ind = 36, num_loci = 13198)

## ------------------------------------------------------------------------
data("mouse_msats")
mouse_msats[1:8, 1:8]

## ------------------------------------------------------------------------
mouse_microsats <- convert_raw(mouse_msats) 
mouse_microsats[1:8, 1:8]

## ---- echo=TRUE, results='hide'------------------------------------------
g2_mouse_microsats <- g2_microsats(mouse_microsats, nperm = 100, nboot = 100, CI = 0.95)
g2_mouse_snps <- g2_snps(mouse_snps, nperm = 100, nboot = 10, 
                         CI = 0.95, parallel = FALSE, ncores = NULL)

## ------------------------------------------------------------------------
g2_mouse_microsats

## ----g2_bootstrapping, echo=TRUE, fig.cap="Distribution of g2 from bootstrapping with confidence interval", fig.width=7, fig.height=4----
par(mfrow=c(1,2))
plot(g2_mouse_microsats, main = "Microsatellites",
     col = "cornflowerblue", cex.axis=0.85)
plot(g2_mouse_snps, main = "SNPs",
     col = "darkgoldenrod1", cex.axis=0.85)

## ---- echo=TRUE, results='hide'------------------------------------------
HHC_mouse_microsats <- HHC(mouse_microsats , niter = 1000)
HHC_mouse_snps <- HHC(mouse_snps, niter = 100)

## ------------------------------------------------------------------------
HHC_mouse_microsats

## ----HHC, echo=TRUE, fig.cap="Distribution of heterozygosity-heterozygosity correlations", fig.width=7, fig.height=4----
par(mfrow=c(1,2))
plot(HHC_mouse_microsats, main = "Microsatellites",
     col = "cornflowerblue", cex.axis=0.85)
plot(HHC_mouse_snps, main = "SNPs",
     col = "darkgoldenrod1", cex.axis=0.85)

## ---- echo=TRUE----------------------------------------------------------
# r^2 between inbreeding and heterozygosity
hf <- r2_hf(genotypes = mouse_microsats, type = "msats")
# r^2 between inbreeding and fitness
Wf <- r2_Wf(genotypes = mouse_microsats, trait = bodyweight, 
            family = gaussian, type = "msats")

## ---- echo=TRUE, results = "hide"----------------------------------------
# r^2 between inbreeding and heterozygosity with bootstrapping
hf <- r2_hf(genotypes = mouse_microsats, nboot = 100, type = "msats", parallel = FALSE)

## ------------------------------------------------------------------------
plot(hf, plottype = "histogram")

## ---- echo=TRUE, results='hide'------------------------------------------
# g2
g2 <- g2_microsats(mouse_microsats)
# calculate sMLH
het <- sMLH(mouse_microsats)
# variance in sMLH
het_var <- var(het)
# Linear model of fitness trait on heterozygosity
mod <- lm(bodyweight ~ het)
# regression slope
beta <- coef(mod)[2]
# r2 between fitness and heterozygosity
Wh <- cor(bodyweight,predict(mod))^2
# r2 between inbreeding and heterozygosity
hf <- r2_hf(genotypes = mouse_microsats, type = "msats")
# r2 between inbreeding and fitness
Wf <- r2_Wf(genotypes = mouse_microsats, trait = bodyweight, 
            family = gaussian, type = "msats")

## ---- echo=FALSE---------------------------------------------------------
library(knitr)
df_msats <- data.frame(g2[["g2"]], het_var, beta, Wh, hf$r2_hf_full, Wf$r2_Wf_full)
row.names(df_msats) <- c("microsats")

kable(df_msats, digits=3, caption="Descriptors of HFCs",
      col.names = c("$\\hat{g}_2$", "$\\hat{\\sigma}^2(h)$", "$\\hat{\\beta}_{Wh}$", "$\\hat{r}^2_{Wh}$", "$\\hat{r}^2_{hf}$","$\\hat{r}^2_{Wf}$"))


## ---- echo=TRUE, results='hide'------------------------------------------
resamp_g2_mouse_microsats <- resample_g2(mouse_microsats, subsets = c(2,4,6,8,10,12), 
                                     nboot = 100, type = "msats")
resamp_g2_mouse_snps <- resample_g2(mouse_snps, subsets = c(100, 500, 1000, 2000, 5000, 
                                    13000), nboot = 10, type = "snps")

## ----resample_g2, echo=TRUE, fig.cap="g2 for different subsets of markers", fig.width=7, fig.height=4----
par(mfrow = c(1, 2))
plot(resamp_g2_mouse_microsats, main = "Microsatellites", col = "cornflowerblue", cex.axis=0.85)
plot(resamp_g2_mouse_snps, main = "SNPs", col = "darkgoldenrod1", cex.axis=0.85)

## ---- echo=TRUE, results='hide'------------------------------------------
r2_hf_mouse_microsats <- r2_hf(mouse_microsats, subsets = c(2,4,6,8,10,12), 
                           nboot = 100, type = "msats")
r2_hf_mouse_snps <- r2_hf(mouse_snps, subsets = c(100, 500, 1000, 2000, 5000, 13000), 
                          nboot = 10, type = "snps")

## ----r2_hf, echo=FALSE, fig.cap="Expected r2 between inbreeding level (f) and heterozygosity", fig.width=7, fig.height=4----
par(mfrow = c(1, 2))
plot(r2_hf_mouse_microsats, main = "Microsatellites", col = "cornflowerblue", cex.axis=0.85)
plot(r2_hf_mouse_snps, main = "SNPs", col = "darkgoldenrod1", cex.axis=0.85)

## ---- results = "hide"---------------------------------------------------
g2_seals <- g2_microsats(mouse_microsats, nperm = 100, 
                         nboot = 100, CI = 0.95)

## ------------------------------------------------------------------------
str(g2_seals)

## ------------------------------------------------------------------------
g2_bootstrap_results <- g2_seals$g2_boot
str(g2_bootstrap_results)

