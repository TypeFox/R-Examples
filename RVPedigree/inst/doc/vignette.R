## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----load.pkg------------------------------------------------------------
library(RVPedigree)

## ----load.pheno----------------------------------------------------------
  pheno.file <- system.file("extdata", "pheno.dat",
                            package="RVPedigree")
  pheno.table <- read.table(pheno.file, header=FALSE)
  dim(pheno.table)
  pheno <- pheno.table[, 1]
  head(pheno)

## ----load.covar----------------------------------------------------------
  covar.file <- system.file("extdata", "covariates.dat",
                            package="RVPedigree")
  covar <- as.matrix(read.table(covar.file, header=FALSE))
  dim(covar)

## ----load.kinmat---------------------------------------------------------
kinship.file <- system.file("extdata", "kinmat.Rdata",
                            package="RVPedigree")
load(kinship.file)
dim(kin1)
rel.mat <- 2 * kin1

## ----relmat.pedigree.ped, eval=FALSE-------------------------------------
#  rel.mat.pedigree.ped <- GetRelMatrix(datatype="pedigree",
#                                       plinkbasefile="mypedfile")

## ----relmat.gen.ped, eval=FALSE------------------------------------------
#  rel.mat.genomic.ped <- GetRelMatrix(datatype="genomic",
#                                      path2Plink="/tools/genepi/Plink_1.90/plink",
#                                      plinkbasefile="mypedfile")

## ----relmat.pedigree.bed, eval=FALSE-------------------------------------
#  rel.mat.pedigree.bed <- GetRelMatrix(datatype="pedigree",
#                                       is.binary=TRUE,
#                                       plinkbasefile="mybedfile")

## ----relmat.gen.bed, eval=FALSE------------------------------------------
#  rel.mat.genomic.bed <- GetRelMatrix(datatype="genomic",
#                                      path2Plink="plink_1.90",
#                                      is.binary=TRUE,
#                                      plinkbasefile="mybedfile")

## ----load.kinmat-GenABEL-------------------------------------------------
genabel.data <- system.file("extdata", "gwaa.data.RData",
                            package="RVPedigree")
library(GenABEL)
load(genabel.data)
nsnps(data1)
nids(data1)

## ----compute.kinmat-GenABEL, eval=FALSE----------------------------------
#  rel.mat.large <- GetRelMatrix("genomic", gwaa.data=data1)

## ----load.map------------------------------------------------------------
map.file <- system.file("extdata", "test.map",
                        package="RVPedigree")
mymap <- readMapFile(map.file)

## ----set.genofile--------------------------------------------------------
geno.file <- system.file("extdata", "test.ped",
                         package="RVPedigree")

## ----ASKAT.region--------------------------------------------------------
askat.result <- ASKAT.region(y=pheno, X=covar, Phi=rel.mat,
                             type="ped",
                             filename=geno.file,
                             chr=20,
                             startpos=0,
                             endpos=10000000,
                             map=mymap,
                             regionname="Gene 1")
print(askat.result)

## ----compute.kinmat-eigenv-----------------------------------------------
relmat.eig    <- eigen(rel.mat)
relmat.eigvec <- relmat.eig$vectors
relmat.eigval <- relmat.eig$values

## ----NASKAT.region-------------------------------------------------------
nASKAT.result <- NormalizedASKAT.region(y=pheno, X=covar, Phi=rel.mat,
                                        type="ped",
                                        filename=geno.file,
                                        chr=20,
                                        startpos=0,
                                        endpos=10000000,
                                        map=mymap,
                                        regionname="Gene 1",
                                        U=relmat.eigvec,
                                        S=relmat.eigval)
print(nASKAT.result)

## ----VCC1.region---------------------------------------------------------
vcc1.result <- VCC1.region(y=pheno, X=covar, Phi=rel.mat,
                           type="ped",
                           filename=geno.file,
                           chr=20,
                           startpos=0,
                           endpos=10000000,
                           map=mymap,
                           regionname="Gene 1")
print(vcc1.result)

## ----VCC2.region---------------------------------------------------------
vcc2.result <- VCC2.region(y=pheno, X=covar, Phi=rel.mat,
                           type="ped",
                           filename=geno.file,
                           chr=20,
                           startpos=0,
                           endpos=10000000,
                           map=mymap,
                           regionname="Gene 1")
print(vcc2.result)

## ----VCC2.region-perm----------------------------------------------------
vcc2.result <- VCC2.region(y=pheno, X=covar, Phi=rel.mat,
                           type="ped",
                           filename=geno.file,
                           chr=20,
                           startpos=0,
                           endpos=10000000,
                           map=mymap,
                           regionname="Gene 1",
                           Nperm=500)
print(vcc2.result)

## ----VCC3.region---------------------------------------------------------
vcc3.result <- VCC3.region(y=pheno, X=covar, Phi=rel.mat,
                           type="ped",
                           filename=geno.file,
                           chr=20,
                           startpos=0,
                           endpos=10000000,
                           map=mymap,
                           regionname="Gene 1")
print(vcc3.result)

## ----VCC3.region-perm----------------------------------------------------
vcc3.result <- VCC3.region(y=pheno, X=covar, Phi=rel.mat,
                           type="ped",
                           filename=geno.file,
                           chr=20,
                           startpos=0,
                           endpos=10000000,
                           map=mymap,
                           regionname="Gene 1",
                           Nperm=200)
print(vcc3.result)

## ----load.genelist-------------------------------------------------------
genes.file <- system.file("extdata", "genes2.lst",
                          package="RVPedigree")
genes <- read.table(genes.file,
                    header=TRUE,
                    stringsAsFactors=FALSE)
genes
colnames(genes) <- c("Name", "Chr", "StartPos", "EndPos")

## ----RVPedigree----------------------------------------------------------
gwresults <- RVPedigree(method="VCC1",
                        y=pheno,
                        X=covar,
                        Phi=rel.mat,
                        regions=genes,
                        filename=geno.file,
                        type='ped',
                        pvalThreshold=1
                        )
gwresults

## ----RVPedigree-VCC3VCC1-------------------------------------------------
gwresults <- RVPedigree(method="VCC1",
                        y=pheno,
                        X=covar,
                        Phi=rel.mat,
                        regions=genes,
                        filename=geno.file,
                        type='ped',
                        Nperm=200,
                        pvalThreshold=0.4,
                        VCC3afterVCC1=TRUE
                        )
gwresults

## ----VCC2.region-parallel------------------------------------------------
vcc2.result <- VCC2.region(y=pheno, X=covar, Phi=rel.mat,
                           type="ped",
                           filename=geno.file,
                           chr=20,
                           startpos=0,
                           endpos=10000000,
                           map=mymap,
                           regionname="Gene 1",
                           Nperm=500,
                           Ncores=2)
print(vcc2.result)

## ----RVPedigree-parallel-------------------------------------------------
gwresults <- RVPedigree(method="VCC1",
                        y=pheno,
                        X=covar,
                        Phi=rel.mat,
                        regions=genes,
                        filename=geno.file,
                        type='ped',
                        Ncores=2,
                        pvalThreshold=1
                        )
gwresults

