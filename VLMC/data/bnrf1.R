## data() sets directory
load(system.file(file.path("xdata", "bnrf1_both.rda"), package = "VLMC"))
bnrf1EB <- as.factor(strsplit(bnrf1EB, "")[[1L]])
bnrf1HV <- as.factor(strsplit(bnrf1HV, "")[[1L]])

## How MM packaged these :
if(FALSE) {
    nuc <- c("a","c","g","t")
    ddir <- "/archives/Data-Collection/ShumwayStoffer-TSA"
    str(bnrf1EB <- paste(nuc[scan(file.path(ddir,"bnrf1ebv.dat"))],collapse=""))
    str(bnrf1HV <- paste(nuc[scan(file.path(ddir,"bnrf1hvs.dat"))],collapse=""))
    setwd("/u/maechler/R/MM/STATISTICS/VLMC/VLMC/data")
    save(bnrf1EB, bnrf1HV, file="bnrf1.rda")
}
