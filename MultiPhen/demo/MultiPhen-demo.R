
library(MultiPhen)
vignette = "../vignettes"

print(mPhen.options("all",descr = TRUE))
print(mPhen.options("misc",descr = TRUE))
options("mPhen.log10p"=FALSE)
pheno.opts = mPhen.options("pheno.input")
pheno = mPhen.readPhenoFiles(paste(vignette,"pheno.txt",sep="/"), opts = pheno.opts)
pheno.opts$mPhen.sep.pheno=" "
pheno.opts$mPhen.numHeaderRows.pheno = 0
pheno.hapmap = mPhen.readPhenoFiles(paste(vignette,"hapmap2.fam",sep="/"), opts = pheno.opts)
print(names(pheno))
print(pheno$limit)
pheno.hapmap$limit$covariates = pheno.hapmap$limit$phenotypes[1]
pheno.hapmap$limit$phenotypes = pheno.hapmap$limit$phenotypes[-1]
opts = mPhen.options("regression")
phenoObject = mPhen.preparePheno(pheno,opts = opts)
phenoObject.hapmap = mPhen.preparePheno(pheno.hapmap,opts = opts)
numPhenos = length(phenoObject$phenN)
numPhenos.hapmap = length(phenoObject.hapmap$phenN)
geno.opts = mPhen.options("geno.input")
opts = mPhen.options("regression")
geno.opts$mPhen.batch = 100
geno.opts$mPhen.format = "GT"
file = paste(vignette,"ALL.chr21.integrated.phase1.v3.20101123.snps.indels.svs.genotypes.extract.impute",sep="/")
file = paste(vignette,"ALL.chr21.integrated.phase1.v3.20101123.snps.indels.svs.genotypes.extract.vcf.gz",sep="/")
file = paste(vignette,"ALL.chr21.integrated.phase1.v3.20101123.snps.indels.svs.genotypes.extract.zip",sep="/")
genoConnection <-mPhen.readGenotypes(file, opts = geno.opts, indiv = rownames(pheno$pheno))
geno <-genoConnection$genoData  
dimg = dim(geno)
file = paste(vignette,"hapmap2",sep="/")
genoConnection.hapmap <-mPhen.readGenotypes(file, opts = geno.opts, indiv = rownames(pheno.hapmap$pheno))
geno.hapmap <-genoConnection.hapmap$genoData  
print(opts)
resultsJoint = mPhen.assoc(geno, phenoObject,  opts = opts)
sigInds = which(resultsJoint$Res[,,numPhenos+1,2]<0.05)
print(resultsJoint$Res[,sigInds,,2])
opts$mPhen.variable.selection=TRUE
opts$mPhen.link.geno="gaussian"
resultsBackwardSelection =mPhen.assoc(geno, phenoObject,  opts=opts)
#print p-values  
print(resultsBackwardSelection$Res[,,,2])
#print betas
print(resultsBackwardSelection$Res[,,,1])

   opts$mPhen.variable.selection=FALSE
   opts$mPhen.JointModel=FALSE 
   opts$mPhen.inverseRegress=FALSE   
   resultsSingle = mPhen.assoc(geno, phenoObject,  opts)

 opts$mPhen.variable.selection=FALSE
 opts$mPhen.JointModel=TRUE
 opts$mPhen.inverseRegress = FALSE    
 resultsJointG = mPhen.assoc(geno, phenoObject,opts)
 opts$mPhen.link.geno="gaussian"
 resultsCCA =  mPhen.cca(geno, phenoObject,opts=opts)
 resDir = "resultsDir" 
 towrite = list(long.txt = TRUE,   wide.txt = TRUE)
 toplot = list(.manh = TRUE,.qq = TRUE,.heatm = TRUE,.fprint = TRUE )
 plotopts = mPhen.options("plot")
 output=mPhen.writeOutput(resultsJoint,output=resDir, geno = geno, towrite = towrite, toplot = toplot, opts = plotopts)
 output=mPhen.writeOutput(resultsBackwardSelection,output=output, geno = geno, towrite = towrite, toplot = toplot,opts = plotopts)


 output=mPhen.writeOutput(resultsSingle,output=output, geno = geno, towrite = towrite, toplot = toplot, opts = plotopts)
 output=mPhen.writeOutput(resultsJointG,output=output, geno = geno, towrite = towrite, toplot = toplot,opts = plotopts)

if(length(dimg)==2){
blockSize = 2 ## size of each block for partitioning correlation
orthogAll = c(0.9,0.5) ## parameters controlling how 'orthogonal' phenotypes are to each other.  
##First entry is orthogonality between blocks, and second is orthoganility within blocks
#parameters must in interval (0,1), with closer to 1 indicating more orthogonal
 covar = mPhen.sampleCovar(numPhenos,blockSize,orthogAll = orthogAll)
  effDir = c(1,1,0,0,0,0)
  total.variance.explained = 0.1  ## total variance explained by all snps.

  snpIndices <- c(1,2)  
  betag = c(1,0.5)
  genoEffect = geno[,snpIndices,drop=F] %*% (betag/(betag %*% betag))

  pheno.sim = mPhen.simulate(genoEffect,dimnames(geno)[[1]], covar,effDir,
			     total.variance.explained, inverse=FALSE, effDirInReverseEigenspace = FALSE)
  mPhen.plotCorrelation(pheno.sim$pheno[,which(pheno.sim$effDir!=0),drop=F],geno[,snpIndices,drop=F],cex=0.5)
phenoObjectSim = mPhen.preparePheno(pheno.sim,opts = opts)
opts = mPhen.options("regression")
resultsJointSim = mPhen.assoc(geno, phenoObjectSim,  opts = opts)
}


