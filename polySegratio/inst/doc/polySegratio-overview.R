### R code from vignette source 'polySegratio-overview.Rnw'

###################################################
### code chunk number 1: polySegratio-overview.Rnw:58-59
###################################################
library(polySegratio)


###################################################
### code chunk number 2: polySegratio-overview.Rnw:62-64
###################################################
op <- options()
options(width=70, digits=4)


###################################################
### code chunk number 3: polySegratio-overview.Rnw:144-150
###################################################
## obtain expected segregation ratios
##  default is one nulliplex parent so type.parents="heterogeneous"

print(unlist(expected.segRatio(2)))
print(unlist(expected.segRatio("Tetraploid")))
print(expected.segRatio("Octa")$ratio)


###################################################
### code chunk number 4: polySegratio-overview.Rnw:173-177
###################################################
## obtain expected segregation ratios with type.parents="homozygous"

print(unlist(expected.segRatio("tetra",type="homoz")))
print(expected.segRatio("Octa",type="homoz")$ratio)


###################################################
### code chunk number 5: polySegratio-overview.Rnw:183-186
###################################################
## obtain expected segregation ratios with odd ploidy level
a <- expected.segRatio(9)
print(a$ratio)


###################################################
### code chunk number 6: polySegratio-overview.Rnw:221-224
###################################################
mark.sim4 <- sim.autoMarkers(4, dose.proportion=c(0.7,0.3), 
                             n.markers=200, n.individuals = 200)
print(mark.sim4)


###################################################
### code chunk number 7: polySegratio-overview.Rnw:229-230
###################################################
plot(mark.sim4)


###################################################
### code chunk number 8: polySegratio-overview.Rnw:257-260
###################################################
miss.sim4 <- addMisclass(mark.sim4, misclass = 0.1)
miss.sim4 <- addMissing(miss.sim4, na.proportion = 0.2)
print(miss.sim4, col=c(1:6))


###################################################
### code chunk number 9: polySegratio-overview.Rnw:265-266
###################################################
plot(miss.sim4, type="all")


###################################################
### code chunk number 10: polySegratio-overview.Rnw:299-308
###################################################
op <- par(mfrow = c(2, 2))	
cmain <- 1.7
plot(sim.autoMarkers(4,c(0.8,0.2)), main="No overdispersion", cex.main=cmain)
plot(sim.autoMarkers(4,c(0.8,0.2), overdisp=TRUE), main="Shape1 = 50", cex.main=cmain)
plot(sim.autoMarkers(4,c(0.8,0.2), overdisp=TRUE, shape1=15), 
     main="Shape1 = 15", cex.main=cmain)
plot(sim.autoMarkers(4,c(0.8,0.2), overdisp=TRUE, shape1=5), 
     main="Shape1 = 5", cex.main=cmain)
par(op)


###################################################
### code chunk number 11: polySegratio-overview.Rnw:347-351
###################################################
## simulated data
a <- sim.autoMarkers(ploidy = 8, c(0.7,0.2,0.09,0.01), n.markers=200, 
                     n.individuals=100)
print(a)


###################################################
### code chunk number 12: polySegratio-overview.Rnw:358-362
###################################################
## summarise chi-squared test vs true
ac <- test.segRatio(a$seg.ratios, ploidy=8, method="chi.squared")
print(ac)
print(addmargins(table(a$true.doses$dosage, ac$dosage, exclude=NULL)))


###################################################
### code chunk number 13: polySegratio-overview.Rnw:375-379
###################################################
## summarise binomial CI vs true
ab <- test.segRatio(a$seg.ratios, ploidy=8, method="bin", alpha=0.01)
print(ab)
print(addmargins(table(a$true.doses$dosage, ab$dosage, exclude=NULL)))


###################################################
### code chunk number 14: polySegratio-overview.Rnw:396-403
###################################################
## imaginary data frame representing ceq marker names read in from
## spreadsheet
x <- data.frame( col1 = c("agc","","","","gct5","","ccc","",""),
                col2 = c(1,3,4,5,1,2,2,4,6))
print(x)
print(makeLabel(x))
print(cbind(x,lab=makeLabel(x, sep=".")))


###################################################
### code chunk number 15: polySegratio-overview.Rnw:412-420
###################################################
p2 <- sim.autoCross(4,
dose.proportion=list(p01=c(0.7,0.3),p10=c(0.7,0.3),
                     p11=c(0.6,0.2,0.2)))
print(p2, row=c(1:5))

ss <- divide.autoMarkers(p2$markers)

print(ss, row=c(1:5))


