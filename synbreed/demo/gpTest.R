#############################################
## Test gpData object and regress
##
##
## author : Hans-Juergen Auinger
## date : 2011 - 11 - 30
##
##############################################

# number of locations
nLoc <- 3
# number of genotypes
nEntry <- 10
# number of replications
nRep <- 2
# phenotypic residuals
pheno <- data.frame(ID = rep(1:nEntry, nRep*nLoc),
                    Trait1 = rnorm(nRep*nLoc*nEntry, sd=3),
                    Trait2 = rnorm(nRep*nLoc*nEntry)*2,
                    Trait3 = rnorm(nRep*nLoc*nEntry, sd=5),
                    loc = rep(1:nLoc, each=nRep*nEntry),
                    wdh = rep(rep(1:nRep, each=nEntry), nLoc))
# individual values
IDvals <- matrix(c(rnorm(nEntry, sd=1),rnorm(nEntry, sd=2),rnorm(nEntry, sd=4)), ncol=3, byrow=FALSE)
# location values
locVals <- matrix(c(rnorm(nLoc, sd=3), rnorm(nLoc, sd=2), rnorm(nLoc, sd=4)), ncol=3, byrow=FALSE)
# replication values
repVals <-matrix(c(rnorm(nLoc*nRep, sd=5), rnorm(nLoc*nRep, sd=5), rnorm(nLoc*nRep, sd=4)), ncol=3, byrow=FALSE)
# sums
Trait1 <- locVals[pheno$loc, 1] + rep(repVals[, 1], each=nEntry) + IDvals[pheno$ID, 1]
Trait2 <- locVals[pheno$loc, 2] + rep(repVals[, 2], each=nEntry) + IDvals[pheno$ID, 2]
Trait3 <- locVals[pheno$loc, 3] + rep(repVals[, 3], each=nEntry) + IDvals[pheno$ID, 3]
pheno$Trait1 <- pheno$Trait1 + Trait1
pheno$Trait2 <- pheno$Trait2 + Trait2
pheno$Trait3 <- pheno$Trait3 + Trait3
pheno$Trait1 <- pheno$Trait1 + ceiling(abs(min(pheno$Trait1))) + 5
pheno$Trait2 <- pheno$Trait2 + ceiling(abs(min(pheno$Trait2))) + 2
pheno$Trait3 <- pheno$Trait3 + ceiling(abs(min(pheno$Trait3))) + 10
# genotypic values
geno <- matrix(sample(c(0:2), 15*nEntry, replace=TRUE), nrow=nEntry, ncol=15)
rownames(geno) <- 1:nEntry
colnames(geno) <- paste("M", 1:15, sep="")
# map infomation
map <- data.frame(row.names = paste("M", 1:15, sep=""), chr = c(rep(1:2, each=6), 2, 2, 2), pos=c((1:6)*30-30, (1:9)*25-25))
# create a gpData object
gpTest <- create.gpData(pheno=pheno, geno=geno, map=map, repeated = c("loc", "wdh"), modCovar = c("loc", "wdh"))
gpTest <- codeGeno(gpTest, label.heter = "1")
summary(gpTest)
# create a realized relationship matrix
kinRel <- kin(gpTest, ret= "realized")
# test gpMod function
testMod1 <- gpMod (gpData=gpTest, kin=kinRel, model="BLUP", trait=1, repl=NULL, markerEffects=TRUE, fixed=~loc, random=~repl)
summary(testMod1)
testMod2 <- gpMod (gpData=gpTest, model="BLUP", trait=1:3, repl=NULL, markerEffects=TRUE, fixed=~loc, random=~repl)
summary(testMod2)
testMod3 <- gpMod (gpData=gpTest, model="BLUP", trait=1:3, repl=NULL, markerEffects=TRUE)
summary(testMod3)

mod1 <- gpMod(gpData=gpTest, kin=kinRel, model="BLUP", trait=2, repl=NULL, markerEffects=FALSE)

mod1a <- gpMod(gpData=gpTest, kin=kinRel, model="BLUP", trait=2, repl=NULL, markerEffects=FALSE)
mod2 <- gpMod(gpData=gpTest, kin=kinRel, model="BLUP", trait=2, repl=NULL, markerEffects=FALSE)
mod3 <- gpMod(gpData=gpTest, kin=kinRel, model="BLUP", trait=1, repl=NULL, markerEffects=FALSE, fixed=~loc, random=~repl)
mod4 <- gpMod(gpData=gpTest, kin=kinRel, model="BLUP", trait=1, repl=NULL, markerEffects=FALSE, fixed=~loc, random=~repl)
