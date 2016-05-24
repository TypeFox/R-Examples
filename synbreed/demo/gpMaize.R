#############################################
## Genomic prediction using the maize data
## using G-BLUP with a comparison to P-BLUP
## 
##
## author : Valentin Wimmer
## date : 2011 - 11 - 25
##
##############################################

# load the maize data
data(maize)

# first overview
summary(maize)

# visualization of the genetic map
plotGenMap(maize)

# recode marker genotypes
maizeC <- codeGeno(maize)

# compute the genomic relationship matrix
# according to Albrecht et al. (2011)
# as all individuals are DH lines being
# fully homozygous and the phenotypes
# were evaluated in a testcross, the 
# formula of Habier et al (2007) must
# be corrected by an additional factor 4
U <- kin(maizeC,ret="realized")
U <- U/2

# compute the pedigree-based kinship
A <- kin(maizeC,ret="kin",DH=maize$covar$DH)

# plot kinship coefficients
plot(A[maize$covar$genotyped,maize$covar$genotyped])
plot(U)
# genomic kinship is different within full-sib families

# extract true breeding values
# (maize is a simulated data set)
tbv <- maize$covar$tbv[maize$covar$genotyped]

# fit genomic prediction models
#modU <- gpMod(maizeC,kin=U,model="BLUP")
#modA <- gpMod(maizeC,kin=A,model="BLUP")

# correlation with true breeding values
#cor(modA$genVal,tbv)
#cor(modU$genVal,tbv)


