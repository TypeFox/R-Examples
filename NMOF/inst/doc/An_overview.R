### R code from vignette source 'An_overview.Rnw'

###################################################
### code chunk number 1: An_overview.Rnw:23-25
###################################################
require("NMOF")
options(continue = " ", digits = 5, max.print = 10, width = 85)


###################################################
### code chunk number 2: An_overview.Rnw:109-110
###################################################
showExample()


###################################################
### code chunk number 3: An_overview.Rnw:118-119 (eval = FALSE)
###################################################
## vignette(package = "NMOF")  ## display vignette titles


###################################################
### code chunk number 4: An_overview.Rnw:121-123
###################################################
x <- vignette(package = "NMOF")
cat(paste(strwrap(x$results[,"Title"], exdent = 2), collapse = "\n"))


###################################################
### code chunk number 5: An_overview.Rnw:144-145 (eval = FALSE)
###################################################
## file.show(system.file("NMOFex/NMOFman.R", package = "NMOF"))


###################################################
### code chunk number 6: An_overview.Rnw:152-159
###################################################
test.rep <- readLines(system.file("unitTests/report.txt", 
                                  package = "NMOF"))
nt <- gsub(".*\\(([0-9]+) checks?\\).*", "\\1",
           test.rep[grep("\\(\\d+ checks?\\)", test.rep)])
cat("Package version  ", gsub("(.*)[.]([0-9]+)$", "\\1-\\2",
                            packageVersion("NMOF")), "\n",
    "Number of tests: ", sum(as.numeric(nt)), sep = "")


###################################################
### code chunk number 7: An_overview.Rnw:184-187 (eval = FALSE)
###################################################
## require("utils")
## bug.report("[NMOF] Unexpected behaviour in function XXX", 
##             maintainer("NMOF"), package = "NMOF")


