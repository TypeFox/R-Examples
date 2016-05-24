### R code from vignette source 'miRtest.Rnw'

###################################################
### code chunk number 1: DatagenerationX
###################################################
#######################################
### Generate random expression data ###
#######################################
# Generate random miRNA expression data of 3 miRNAs
# with 8 replicates
set.seed(1)
X = rnorm(24);
dim(X) = c(3,8);
rownames(X) = 1:3;


###################################################
### code chunk number 2: DatagenerationY
###################################################
# Generate random mRNA expression data with 20 mRNAs
# and 10 replicates
Y = rnorm(200);
dim(Y) = c(20,10);
rownames(Y) = 1:20;


###################################################
### code chunk number 3: groups
###################################################
# Let's assume that we want to compare 2 miRNA groups, each of 4 replicates:
group.miRNA = factor(c(1,1,1,1,2,2,2,2));
# ... and that the corresponding mRNA experiments had 5 replicates in each group
group.mRNA = factor(c(1,1,1,1,1,2,2,2,2,2));


###################################################
### code chunk number 4: allocation
###################################################
####################
### Perform Test ###
####################
library(miRtest)
#Let miRNA 1 attack mRNAs 1 to 9 and miRNA 2 attack mRNAs 10 to 17.
# mRNAs 18 to 20 are not attacked. miRNA 3 has no gene set.
miR = c(rep(1,9),c(rep(2,8)));
mRNAs = 1:17;
A = data.frame(mRNAs,miR); # Note that the miRNAs MUST be in the second column!
A


###################################################
### code chunk number 5: analysis
###################################################
set.seed(1)
P = miR.test(X,Y,A,group.miRNA,group.mRNA)
P




###################################################
### code chunk number 6: otherTests
###################################################
#####################################################
### For a faster result: use other gene set tests ###
#####################################################
# Wilcoxon two-sample test is recommended for fast results
# Note that results may vary depending on how much genes correlate

P.gsWilcox = miR.test(X,Y,A,group.miRNA,group.mRNA,gene.set.tests="W")
P.gsWilcox


###################################################
### code chunk number 7: otherA
###################################################
############################################
### We can use an allocation matrix as A ###
############################################
A = generate.A(A,X=X,Y=Y,verbose=FALSE);
A


###################################################
### code chunk number 8: otherAtest
###################################################
# Now we can test as before
set.seed(1)
P = miR.test(X,Y,A,group.miRNA,group.mRNA,allocation.matrix=TRUE)
P




###################################################
### code chunk number 9: otherDesigns
###################################################
#####################
### Other Designs ###
#####################

# Some more complicated designs are implemented, check the vignette "miRtest" for details.
group.miRNA = 1:8
group.mRNA = 1:10
covariable.miRNA = factor(c(1,2,3,4,1,2,3,4))    ### A covariable in miRNAs.
covariable.mRNA = factor(c(1,2,3,4,5,1,2,3,4,5)) ### A covariable in mRNAs.

library(limma)
design.miRNA = model.matrix(~group.miRNA + covariable.miRNA)
design.mRNA =  model.matrix(~group.mRNA + covariable.mRNA)



###################################################
### code chunk number 10: otherDesignsMirTest
###################################################
P = miR.test(X,Y,A,design.miRNA=design.miRNA,design.mRNA=design.mRNA,allocation.matrix=TRUE)
P


###################################################
### code chunk number 11: example
###################################################
library(miRtest)
data("X",package="miRtest")  ### this is the miRNA expression data of the ten most interesting miRNAs
data("Y",package="miRtest")  ### this is the mRNA expression data of the genes targeted by these miRNAs
data("A",package="miRtest")  ### allocation data from TargetScan


###################################################
### code chunk number 12: example2
###################################################
group.mRNA = factor(c(1,1,1,1,2,2,2,2));
group.miRNA = factor(c(1,1,1,2,2,2));


###################################################
### code chunk number 13: example3
###################################################
P = miR.test(X=X,Y=Y,A=A,group.mRNA=group.mRNA,group.miRNA=group.miRNA,
adjust="BH",gene.set.tests="all",verbose=TRUE)


