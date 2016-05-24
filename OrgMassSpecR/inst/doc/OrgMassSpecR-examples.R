
## ----Setup, echo=FALSE---------------------------------------------------
library(knitr)
opts_chunk$set(comment=NA)


## ----Load----------------------------------------------------------------
library(OrgMassSpecR)


## ----Example1------------------------------------------------------------
MonoisotopicMass(formula = list(C=14, H=8, Cl=4))


## ----Example2------------------------------------------------------------
MonoisotopicMass(formula = list(C=14, H=8, Cl=3))
MonoisotopicMass(formula = list(C=14, H=8, Cl=2))


## ----Example3------------------------------------------------------------
MonoisotopicMass(formula = list(C=2, H=8, Cl=4, x = 12), 
                 isotopes = list(x = 13.0033548378))


## ----Example4------------------------------------------------------------
dde.dist <- IsotopicDistribution(formula = list(C=14, H=8, Cl=4))
dde.dist


## ----Example5, fig.width=4, fig.height=4, dpi=300, out.width="400px", out.height="400px"----
# plot
library(lattice)
print(xyplot(percent ~ mz,
  data = dde.dist,
  type = "h",
  xlab = "m/z",
  ylab = "intensity (%)",
  main = "Isotopic Distribution, DDE")
)


## ----Example6------------------------------------------------------------
hsa <- Digest(example.sequence)
head(hsa)


## ----Example7------------------------------------------------------------
hsa.sub <- subset(hsa, nchar(hsa$peptide) >= 5 & nchar(hsa$peptide) <= 12)
head(hsa.sub)


## ----Example8------------------------------------------------------------
transitions <- FragmentPeptide(c("YLYEIAR", "AEFAEVSK"))
head(transitions)


## ----Example9------------------------------------------------------------
c13.labeled <- FragmentPeptide("YLYEIAr", custom = list(code = "r", 
  mass = MonoisotopicMass(formula = list(C=6, H=12, N=4, O=1), 
                          isotopes = list(C=13.0033548378))))
head(c13.labeled)


## ----Example10-----------------------------------------------------------
n15.labeled <- FragmentPeptide("YLYEIAR", N15 = TRUE)
head(n15.labeled)


## ----Example11, fig.width=4, fig.height=4, dpi=300, out.width="400px", out.height="400px"----
theoretical.dist <- IsotopicDistributionN("YEVQGEVFTKPQLWP", incorp = 0.99)
print(xyplot(percent ~ mz,
  data = theoretical.dist,
  type = "h",
  xlab = "m/z",
  ylab = "intensity (%)",
  main = "Theoretical Isotopic Distribution,\n YEVQGEVFTKPQLWP, 99% 15N")
)


## ----Example12, fig.width=4, fig.height=4, dpi=300, out.width="400px", out.height="400px"----
example.spectrum.labeled$percent <- with(example.spectrum.labeled, 
  intensity / max(intensity) * 100)
print(xyplot(percent ~ mz,
  data = example.spectrum.labeled,
  type = "l",
  xlab = "m/z",
  ylab = "intensity (%)",
  main = "Measured Isotopic Distribution,\n YEVQGEVFTKPQLWP")
)


