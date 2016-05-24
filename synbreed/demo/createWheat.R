#############################################
## Create a genomic prediction Data object
## from the wheat data in the BLR package
##
##
## author : Valentin Wimmer
## date : 2011 - 11 - 30
##
##############################################

data(wheat)
# X = genotypes
# Y = phenotypes

# assign names
rownames(X) <- colnames(A)
colnames(Y) <- paste("Env",1:4,sep="")
rownames(Y) <- colnames(A)

# create a gpData object
gpWheat <- create.gpData(pheno=Y,geno=X)
gpWheat <- codeGeno(gpWheat)

# predictive ability using Bayesian Lasso
# use prior values from Crossa et al. (2010)
priorCrossa <- list(varE=list(df=4,S=1),lambda = list(shape=0.6,rate=1e-4,value=20,type='random'))

# model M-BL for GY-E1 
modMBL <- gpMod(gpWheat,trait=1,model="BL",prior=priorCrossa,nIter=15000,burnIn=5000)

# model PM-BL for GY-E1
modPMBL <- gpMod(gpWheat,trait=1,model="BL",prior=priorCrossa,nIter=15000,burnIn=5000,kin=A)


# extract predicted genetic values and plot versus phenotypic values
plot(predict(modPMBL),gpWheat$pheno[,1,])

cv <- crossVal(gpWheat,trait=1,VC.est="BL",prior=priorCrossa,k=10,Rep=1,Seed=123)
summary(cv)
#Object of class 'cvData'
#
# 10 -fold cross validation with 1 replication(s)
#     Sampling:                 random
#     Variance components:      reestimated with BL
#     Number of random effects: 1279
#     Number of individuals:    599 -- 599
#     Size of the TS:           59 -- 60
#
#Results:
#                      Min         Mean +- pooled SE       Max
# Predictive ability:  0.4035      0.5341 +- NA            0.6614
# Rank correlation:    0.3519      0.4765 +- NA            0.6283
# Mean squared error:  0.492       0.727 +- NA     0.901
# Bias:                0.8252      1.0393 +- NA            1.2426
# 10% best predicted:  0.58        0.58 +- NA      0.58
#
#Seed start:  123
#Seed replications:
#[1] 28758

