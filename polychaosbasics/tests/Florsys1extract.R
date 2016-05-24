###########################################################
# Computation of Sensitivity Indexes (SI)
# for the dataset FLORSYS1.txt 
# via Polynomial Chaos Expansion (PCE)
###########################################################


# Load of necessary functions
library("polychaosbasics")


degree<-4
#############################################
# Read data
#############################################
nomlhs<-'FLORSYS1.txt'
lhdnat<-read.table(system.file("extdata", nomlhs, package="polychaosbasics"))

#lhdnat<-read.table(paste("../inst/extdata/", nomlhs, sep=""))


#############################################
# Extract some data for reducing execution time
#############################################
 nlhs <- 200 # nlhs <- nrow(lhd0)
nvx <- 3
Y<-lhdnat[1:nlhs,ncol(lhdnat)]
lhdnat<-lhdnat[1:nlhs, 1:nvx, drop=FALSE]

#############################################
# Build Legendre polynomial
#############################################
pce <- polyLeg(lhdnat, Y, degree)

# -----------------------------------------
# Default display
# -----------------------------------------
print(pce)

# -----------------------------------------
# All the components in the returned object
# -----------------------------------------
getNames(pce)
# -----------------------------------------
# PCEdesign object
# -----------------------------------------
print(pce@design, all=TRUE)

#############################################
# PCESI calculation
#############################################
retour<- PCESI(pce)
print(retour, all=TRUE)
# -----------------------------------------
# All the components in the returned object
# -----------------------------------------
getNames(retour)


#############################################
# Plot Y against Y.hat
#############################################
jpeg(file="FLORSYS1.jpg")
 plot(pce[, "Y"], retour@y.hat, xlab="output", ylab="predicted",
      main="FLORSYS1 test")
 dev.off()

