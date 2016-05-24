### R code from vignette source 'kinfit.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(prompt = "R> ")
options(SweaveHooks = list(
  cex = function() par(cex.lab = 1.3, cex.axis = 1.3)))


###################################################
### code chunk number 2: FOCUS_2006_C_data
###################################################
library("kinfit")
data("FOCUS_2006_C", package = "kinfit")
print(FOCUS_2006_C)


###################################################
### code chunk number 3: data_format
###################################################
kindata_example <- data.frame(
  t = c(0, 1, 3, 7, 14, 28, 63, 91, 119),
  parent = c(85.1, 57.9, 29.9, 14.6, 9.7, 6.6, 4, 3.9, 0.6)
)


###################################################
### code chunk number 4: FOCUS_2006_C_fits
###################################################
kinfits.C <- kinfit(FOCUS_2006_C, kinmodels = c("SFO", "FOMC", "DFOP", "HS"))


###################################################
### code chunk number 5: FOCUS_2006_C_results
###################################################
kinresults(kinfits.C)


###################################################
### code chunk number 6: FOCUS_2006_C_object
###################################################
kinobject.C <- kinobject <- list(
        parent = "Compound XY",
        type = "Degradation in the environment",
        system = "System 1",    
        source = "Synthetic example data from FOCUS kinetics",
        data = FOCUS_2006_C,
        fits = kinfits.C,
        results = kinresults(kinfits.C))


###################################################
### code chunk number 7: FOCUS_2006_C_report
###################################################
kinreport(kinobject.C)


###################################################
### code chunk number 8: FOCUS_2006_C_figure
###################################################
getOption("SweaveHooks")[["cex"]]()
kinfits.C <- kinfit(FOCUS_2006_C, kinmodels = c("SFO", "FOMC", "DFOP", "HS"))
kinplot(kinobject.C)
title("FOCUS dataset C")


###################################################
### code chunk number 9: FOCUS_2006_C_res
###################################################
getOption("SweaveHooks")[["cex"]]()
kinresplot(kinobject.C, "SFO")


###################################################
### code chunk number 10: kinfits_A_to_F
###################################################
datasets <- LETTERS[1:4]
data(list=paste("FOCUS_2006_", datasets, sep=""), package = "kinfit")
kinobjects <- list()
for (dataset in datasets)
{
  kinobjects[[dataset]] <- list()
  kinobjects[[dataset]]$data <- get(paste("FOCUS_2006_", dataset, sep=""))
  kinobjects[[dataset]]$fits <- 
    kinfit(kinobjects[[dataset]]$data,
    kinmodels = c("SFO", "FOMC", "DFOP", "HS"))
  kinobjects[[dataset]]$results <- 
    kinresults(kinobjects[[dataset]]$fits)
}

data(FOCUS_2006_F, package = "kinfit")
# Set the initial concentration in the sediment to zero
FOCUS_2006_F[1, "parent.sediment"] <- 0
# Calculate total system values for the parent compound
FOCUS_2006_F = transform(FOCUS_2006_F, 
  parent.system = parent.water + parent.sediment)

subsets <- c("system", "water")
for (subset in subsets)
{
  index <- paste("F", subset, sep=" ")
  kinobjects[[index]] <- list()
  kinobjects[[index]]$data <- data.frame(
    t = FOCUS_2006_F$t,
    parent = FOCUS_2006_F[[paste("parent", subset, sep=".")]])
  kinobjects[[index]]$fits <- 
    kinfit(kinobjects[[index]]$data,
    kinmodels = c("SFO", "FOMC", "DFOP", "HS"))
  kinobjects[[index]]$results <- 
    kinresults(kinobjects[[index]]$fits)
}


###################################################
### code chunk number 11: SFO
###################################################
data("FOCUS_2006_SFO_ref_A_to_F", package = "kinfit")
kinmodel = "SFO"
refs <- list()
for (kinobjectname in names(kinobjects))
{
  ref <- subset(FOCUS_2006_SFO_ref_A_to_F, dataset == kinobjectname)
  ref <- ref[-6]
  texfile <- paste("FOCUS_2006_", kinmodel, "_", 
    gsub(" ", "_", kinobjectname), "_ref.tex", sep="")
  write.table(format(ref, nsmall=2), 
    file = texfile,
    sep=" & ", quote=FALSE, 
    row.names=FALSE, col.names=FALSE, eol = " \\\\ \n")
  refs[[kinobjectname]] <- ref
}


###################################################
### code chunk number 12: FOMC
###################################################
data("FOCUS_2006_FOMC_ref_A_to_F", package = "kinfit")
kinmodel = "FOMC"
refs <- list()
for (kinobjectname in names(kinobjects)[c(1:3, 5:6)])
{
  ref <- subset(FOCUS_2006_FOMC_ref_A_to_F, dataset == kinobjectname)
  ref$package <- gsub("#", "$^a$", ref$package)
  ref <- ref[-7]
  texfile <- paste("FOCUS_2006_", kinmodel, "_", 
    gsub(" ", "_", kinobjectname), "_ref.tex", sep="")
  write.table(format(ref, nsmall=2), 
    file = texfile,
    sep=" & ", quote=FALSE, 
    row.names=FALSE, col.names=FALSE, eol = " \\\\ \n")
  refs[[kinobjectname]] <- ref
}


###################################################
### code chunk number 13: DFOP
###################################################
data("FOCUS_2006_DFOP_ref_A_to_B", package = "kinfit")
kinmodel = "DFOP"
refs <- list()
for (kinobjectname in names(kinobjects)[1:2])
{
  ref <- subset(FOCUS_2006_DFOP_ref_A_to_B, dataset == kinobjectname)
  ref <- ref[-8]
  texfile <- paste("FOCUS_2006_", kinmodel, "_", 
    gsub(" ", "_", kinobjectname), "_ref.tex", sep="")
  write.table(format(ref, nsmall=2), 
    file = texfile,
    sep=" & ", quote=FALSE, 
    row.names=FALSE, col.names=FALSE, eol = " \\\\ \n")
  refs[[kinobjectname]] <- ref
}


###################################################
### code chunk number 14: HS
###################################################
data("FOCUS_2006_HS_ref_A_to_F", package = "kinfit")
kinmodel = "HS"
refs <- list()
for (kinobjectname in names(kinobjects)[c(1:3, 5:6)])
{
  ref <- subset(FOCUS_2006_HS_ref_A_to_F, dataset == kinobjectname)
  ref$package <- gsub("\\*", "$^a$", ref$package)
  ref <- ref[-8]
  texfile <- paste("FOCUS_2006_", kinmodel, "_", 
    gsub(" ", "_", kinobjectname), "_ref.tex", sep="")
  write.table(format(ref, nsmall=2), 
    file = texfile,
    sep=" & ", quote=FALSE, 
    row.names=FALSE, col.names=FALSE, eol = " \\\\ \n")
  refs[[kinobjectname]] <- ref
}


###################################################
### code chunk number 15: HS2
###################################################
kinobjects$A$fits.2 <- kinfit(kinobjects$A$data,
  kinmodels = c("HS"),
  start.HS = list(parent.0 = 100, tb = 5, k1 = 0.017, k2 = 0.05))
kinobjects$A$results.2 <- kinresults(kinobjects$A$fits.2)

kinobjects$A$fits.3 <- kinfit(kinobjects$A$data,
  kinmodels = c("HS"),
  start.HS = list(parent.0 = 100, tb = 11, k1 = 0.017, k2 = 0.05))
kinobjects$A$results.3 <- kinresults(kinobjects$A$fits.3)


###################################################
### code chunk number 16: KinGUI_write
###################################################
chi2.SFO.kinfit <- chi2.FOMC.kinfit <- array(dim = length(kinobjects), 
  dimnames = list(names(kinobjects)))
chi2.DFOP.kinfit <- chi2.HS.kinfit <- array(dim = length(kinobjects),
  dimnames = list(names(kinobjects)))
for (kinobjectname in names(kinobjects))
{
  outname <- paste("KinGUI/", gsub(" ", "_", kinobjectname), "_KinGUI.txt", 
    sep="")
  kinwrite.KinGUI(kinobjects[[kinobjectname]], outname)
  chi2.SFO.kinfit[[kinobjectname]] <- 
    kinobjects[[kinobjectname]]$results$stats[["SFO", "err.min"]]
  chi2.FOMC.kinfit[[kinobjectname]] <- 
    ifelse(class(kinobjects[[kinobjectname]]$fits$FOMC) == "try-error",
      NA, kinobjects[[kinobjectname]]$results$stats[["FOMC", "err.min"]])
  chi2.DFOP.kinfit[[kinobjectname]] <- 
    ifelse(class(kinobjects[[kinobjectname]]$fits$DFOP) == "try-error",
      NA, kinobjects[[kinobjectname]]$results$stats[["DFOP", "err.min"]])
  chi2.HS.kinfit[[kinobjectname]] <- 
    ifelse(class(kinobjects[[kinobjectname]]$fits$HS) == "try-error",
      NA, kinobjects[[kinobjectname]]$results$stats[["HS", "err.min"]])
}

chi2.SFO.KinGUI <- c(8.3852, 4.4562, 15.8456, 6.4539, 12.5386, 10.8069)
chi2.FOMC.KinGUI <- c(9.3116, 4.6641, 6.6574, 6.8080, 13.4533, 11.6682)
chi2.DFOP.KinGUI <- c(9.6600, 4.9562, 2.6613, 7.2751, 14.1524, 12.1821)
chi2.HS.KinGUI <- c(4.1106, 4.4535, 4.6963, 5.8196, 3.2178, 1.6558)
names(chi2.SFO.KinGUI) <- names(chi2.FOMC.KinGUI) <- names(kinobjects)
names(chi2.DFOP.KinGUI) <- names(chi2.HS.KinGUI) <- names(kinobjects)

chi2 <- data.frame(
  SFO.KinGUI = chi2.SFO.KinGUI,
  SFO.kinfit = round(100 * chi2.SFO.kinfit, 4),
  FOMC.KinGUI = chi2.FOMC.KinGUI,
  FOMC.kinfit = round(100 * chi2.FOMC.kinfit, 4),
  DFOP.KinGUI = chi2.DFOP.KinGUI,
  DFOP.kinfit = round(100 * chi2.DFOP.kinfit, 4),
  HS.KinGUI = chi2.HS.KinGUI,
  HS.kinfit = round(100 * chi2.HS.kinfit, 4)
)
write.table(chi2,
  file = "chi2_comparison.tex",
    sep=" & ", quote=FALSE, na="",
    row.names=TRUE, col.names=FALSE, eol = " \\\\ \n")


###################################################
### code chunk number 17: KinGUI_write
###################################################
R2.SFO.kinfit <- R2.FOMC.kinfit <- array(dim = 4, 
  dimnames = list(names(kinobjects[1:4])))
R2.DFOP.kinfit <- R2.HS.kinfit <- array(dim = 4,
  dimnames = list(names(kinobjects[1:4])))
for (kinobjectname in names(kinobjects[1:4]))
{
  R2.SFO.kinfit[[kinobjectname]] <- 
    kinobjects[[kinobjectname]]$results$stats[["SFO", "R2"]]
  R2.FOMC.kinfit[[kinobjectname]] <- 
    ifelse(class(kinobjects[[kinobjectname]]$fits$FOMC) == "try-error",
      NA, kinobjects[[kinobjectname]]$results$stats[["FOMC", "R2"]])
  R2.DFOP.kinfit[[kinobjectname]] <- 
    ifelse(class(kinobjects[[kinobjectname]]$fits$DFOP) == "try-error",
      NA, kinobjects[[kinobjectname]]$results$stats[["DFOP", "R2"]])
  R2.HS.kinfit[[kinobjectname]] <- 
    ifelse(class(kinobjects[[kinobjectname]]$fits$HS) == "try-error",
      NA, kinobjects[[kinobjectname]]$results$stats[["HS", "R2"]])
}

EF.SFO.KinGUI <- c(0.9845, 0.9971, 0.9714, 0.9919)
EF.FOMC.KinGUI <- c(0.9831, 0.9973, 0.9955, 0.9920)
EF.HS.KinGUI <- c(0.9972, 0.9972, 0.9980, 0.9945)
EF.DFOP.KinGUI <- c(0.9845, 0.9973, 0.9994, 0.9919)
names(EF.SFO.KinGUI) <- names(EF.FOMC.KinGUI) <- names(kinobjects[1:4])
names(EF.DFOP.KinGUI) <- names(EF.HS.KinGUI) <- names(kinobjects[1:4])

R2 <- data.frame(
  SFO.KinGUI = EF.SFO.KinGUI,
  SFO.kinfit = round(R2.SFO.kinfit, 4),
  FOMC.KinGUI = EF.FOMC.KinGUI,
  FOMC.kinfit = round(R2.FOMC.kinfit, 4),
  DFOP.KinGUI = EF.DFOP.KinGUI,
  DFOP.kinfit = round(R2.DFOP.kinfit, 4),
  HS.KinGUI = EF.HS.KinGUI,
  HS.kinfit = round(R2.HS.kinfit, 4)
)
write.table(R2,
  file = "R2_comparison.tex",
    sep=" & ", quote=FALSE, na="",
    row.names=TRUE, col.names=FALSE, eol = " \\\\ \n")


