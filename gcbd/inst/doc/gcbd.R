### R code from vignette source 'gcbd.Rnw'

###################################################
### code chunk number 1: prelim
###################################################
options(width=50)
library(gcbd)
gcbd.version <- packageDescription("gcbd")$Version
gcbd.date <- packageDescription("gcbd")$Date
now.date <- strftime(Sys.Date(), "%B %d, %Y")
# create figures/ if not present
if ( ! (file.exists("figures") && file.info("figures")$isdir) ) dir.create("figures")


###################################################
### code chunk number 2: data
###################################################
i7 <- getBenchmarkData("i7_920")
xeon <- getBenchmarkData("xeon_X5570")
D <- subset(i7[,-c(1:2,5)], type=="matmult")


###################################################
### code chunk number 3: Lattice
###################################################
figure_Lattice(titles=FALSE)


###################################################
### code chunk number 4: LogLogLattice
###################################################
figure_LogLogLattice(titles=FALSE)


###################################################
### code chunk number 5: table
###################################################
C <- loglogAnalysis()[["coefs"]]
C[,-c(1,2)] <- format(round(C[,-c(1,2)], digits=2))


###################################################
### code chunk number 6: loglogslopes
###################################################
figure_LogLogSlopes()


