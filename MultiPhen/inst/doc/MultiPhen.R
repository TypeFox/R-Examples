### R code from vignette source 'MultiPhen.Rnw'

###################################################
### code chunk number 1: MultiPhen.Rnw:53-54
###################################################
R.version.string


###################################################
### code chunk number 2: MultiPhen.Rnw:58-59 (eval = FALSE)
###################################################
## install.packages("MultiPhen", dependecies = TRUE)


###################################################
### code chunk number 3: MultiPhen.Rnw:63-64
###################################################
library(MultiPhen)


###################################################
### code chunk number 4: MultiPhen.Rnw:74-79
###################################################
library(MultiPhen)
options(width = 60) # for printing this vignette is a reasonably acceptable way
print(mPhen.options("all",descr = TRUE))
print(mPhen.options("misc",descr = TRUE))
options("mPhen.log10p"=FALSE)


###################################################
### code chunk number 5: MultiPhen.Rnw:86-91
###################################################
pheno.opts = mPhen.options("pheno.input")
pheno = mPhen.readPhenoFiles("pheno.txt", opts = pheno.opts)
pheno.opts$mPhen.sep.pheno=" "
pheno.opts$mPhen.numHeaderRows.pheno = 0
pheno.hapmap = mPhen.readPhenoFiles("hapmap2.fam", opts = pheno.opts)


###################################################
### code chunk number 6: MultiPhen.Rnw:96-100
###################################################
print(names(pheno))
print(pheno$limit)
pheno.hapmap$limit$covariates = pheno.hapmap$limit$phenotypes[1]
pheno.hapmap$limit$phenotypes = pheno.hapmap$limit$phenotypes[-1]


###################################################
### code chunk number 7: MultiPhen.Rnw:107-112
###################################################
opts = mPhen.options("regression")
phenoObject = mPhen.preparePheno(pheno,opts = opts)
phenoObject.hapmap = mPhen.preparePheno(pheno.hapmap,opts = opts)
numPhenos = length(phenoObject$phenN)
numPhenos.hapmap = length(phenoObject.hapmap$phenN)


###################################################
### code chunk number 8: MultiPhen.Rnw:118-122
###################################################
geno.opts = mPhen.options("geno.input")
opts = mPhen.options("regression")
geno.opts$mPhen.batch = 100
geno.opts$mPhen.format = "GT"


###################################################
### code chunk number 9: MultiPhen.Rnw:129-138
###################################################
file = "ALL.chr21.integrated.phase1.v3.20101123.snps.indels.svs.genotypes.extract.zip"
genoConnection <-mPhen.readGenotypes(file, opts = geno.opts, 
                                     indiv = rownames(pheno$pheno))
geno <-genoConnection$genoData  
dimg = dim(geno)
file = "hapmap2"
genoConnection.hapmap <-mPhen.readGenotypes("hapmap2", opts = geno.opts, 
                                            indiv = rownames(pheno.hapmap$pheno))
geno.hapmap <-genoConnection.hapmap$genoData  


###################################################
### code chunk number 10: MultiPhen.Rnw:146-148
###################################################
print(opts)
resultsJoint = mPhen.assoc(geno, phenoObject,  opts = opts)


###################################################
### code chunk number 11: MultiPhen.Rnw:154-156
###################################################
sigInds = which(resultsJoint$Res[,,numPhenos+1,2]<0.05)
print(resultsJoint$Res[,sigInds,,2])


###################################################
### code chunk number 12: MultiPhen.Rnw:161-169
###################################################
opts$mPhen.variable.selection=TRUE
opts$mPhen.link.geno="gaussian"
resultsBackwardSelection =mPhen.assoc(geno, 
                                      phenoObject,  opts=opts)
#print p-values  
print(resultsBackwardSelection$Res[,,,2])
#print betas
print(resultsBackwardSelection$Res[,,,1])


###################################################
### code chunk number 13: MultiPhen.Rnw:177-181
###################################################
   opts$mPhen.variable.selection=FALSE
   opts$mPhen.JointModel=FALSE 
   opts$mPhen.inverseRegress=FALSE   
   resultsSingle = mPhen.assoc(geno, phenoObject,  opts)


###################################################
### code chunk number 14: MultiPhen.Rnw:186-190
###################################################
  opts$mPhen.variable.selection=FALSE
  opts$mPhen.JointModel=TRUE
  opts$mPhen.inverseRegress = FALSE    
  resultsJointG = mPhen.assoc(geno, phenoObject,opts)


###################################################
### code chunk number 15: MultiPhen.Rnw:198-201
###################################################
opts$mPhen.link.geno="gaussian"
resultsCCA =  mPhen.cca(geno, 
                        phenoObject,opts=opts)


###################################################
### code chunk number 16: MultiPhen.Rnw:209-217
###################################################
 resDir = "resultsDir" 
 towrite = list(long.txt = TRUE,   wide.txt = TRUE)
 toplot = list(.manh = TRUE,.qq = TRUE,.heatm = TRUE,
               .fprint = TRUE)
 plotopts = mPhen.options("plot")
 output=mPhen.writeOutput(resultsJoint,output=resDir, 
                          geno = geno, towrite = towrite, 
                          toplot = toplot, opts = plotopts)


###################################################
### code chunk number 17: MultiPhen.Rnw:222-231
###################################################
 output=mPhen.writeOutput(resultsBackwardSelection,output=output, 
                          geno = geno, towrite = towrite, 
                          toplot = toplot,     opts = plotopts)
 output=mPhen.writeOutput(resultsSingle,output=resDir, geno = geno, 
                          towrite = towrite, toplot = toplot, 
                          opts = plotopts)
 output=mPhen.writeOutput(resultsJointG,output=output,
                          geno = geno, towrite = towrite,
                          toplot = toplot,     opts = plotopts)


###################################################
### code chunk number 18: MultiPhen.Rnw:241-248
###################################################
blockSize = 2 ## size of each block for partitioning correlation
orthogAll = c(0.9,0.5) ## parameters controlling 
## how 'orthogonal' phenotypes are to each other.  
## First entry is orthogonality between blocks, and 
## second is orthoganility within blocks
#parameters must in interval (0,1), with closer to 1 indicating more orthogonal
 covar = mPhen.sampleCovar(numPhenos,blockSize,orthogAll = orthogAll)


###################################################
### code chunk number 19: MultiPhen.Rnw:253-265
###################################################
  effDir = c(1,1,0,0,0,0)
  total.variance.explained = 0.1  ## total variance explained by all snps.

  snpIndices <- c(1,2)  
  betag = c(1,0.5)
  genoEffect = geno[,snpIndices,drop=F] %*% (betag/(betag %*% betag))

  pheno.sim = mPhen.simulate(genoEffect,dimnames(geno)[[1]], covar,effDir,
			     total.variance.explained, inverse=FALSE, 
           effDirInReverseEigenspace = FALSE)
  mPhen.plotCorrelation(pheno.sim$pheno[,which(pheno.sim$effDir!=0),
                        drop=F],geno[,snpIndices,drop=F],cex=0.5)


###################################################
### code chunk number 20: MultiPhen.Rnw:270-273
###################################################
phenoObjectSim = mPhen.preparePheno(pheno.sim,opts = opts)
opts = mPhen.options("regression")
resultsJointSim = mPhen.assoc(geno, phenoObjectSim,  opts = opts)


