### R code from vignette source 'acnr.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: acnr.Rnw:22-23
###################################################
library(acnr)


###################################################
### code chunk number 2: acnr.Rnw:29-48
###################################################
tf <- .5
dataSets <- c("GSE29172", "GSE11976")

regList <- lapply(dataSets, FUN=function(ds) {
    regDat <- loadCnRegionData(dataSet=ds, tumorFraction=tf) 
    regs <- regDat[["region"]]
})
names(regList) <- dataSets

allregs <- unique(unlist(regList))
tab <- sapply(allregs, FUN=function(reg) {
    sapply(regList, FUN=function(rr) sum(rr==reg))
})
cap <- paste("Size of annotated copy-number regions for each of the", 
            length(dataSets), "data sets.")
if (require(xtable)) {
    print(xtable(tab, caption=cap, label="tab:regData",
                table.placement="!h", caption.placement="bottom"))
}


###################################################
### code chunk number 3: acnr.Rnw:52-53
###################################################
sessionInfo()


