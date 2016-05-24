## This demo will give a typical workflow
##
##
library(seqminer)

## ---------- Step 1 ----------
## Annotate gene
## We prepare annotation parameters:
##   referecne specifies reference genome
##   geneFile is a refFlat format gene definition file
param <- list(reference = system.file("tabanno/test.fa", package = "seqminer"),
              geneFile = system.file("tabanno/test.gene.txt", package = "seqminer"))
param <- makeAnnotationParameter(param)
pause <- readline("Hit ENTER to continue...")


## Then we use an example VCF file and let it output to your current directory.
## The output file is named as "out.vcf.gz"
input <- system.file("tabanno/input.demo.vcf", package = "seqminer")
output <- paste0(getwd(), "/", "out.vcf.gz")
annotateVcf (input, output, param)
pause <- readline("Hit ENTER to continue...")

## ---------- Step 1 Advance ----------
## Integrate external bioinformatics DB
## Here we use a region-based file as example (BED format).
## We can combine the gene annotation and DB integration in one step:
## Note,
##  1. we use 'bed=' option to specify a region-based database;
##  2. we use 'tabix=' option to specify a tabix-based resource file;
##  3. we use 'indexOutput=TRUE' to index outputted VCF file
setwd(system.file("tabanno", package = "seqminer"))
param  <- list(reference = "test.fa",
               geneFile = "test.gene.txt",
               bed = "REGION=test.bed",
               tabix = "test.dbNSFP.gz(SIFT=9,PolyPhen=10)",
               indexOutput = TRUE)
param  <- makeAnnotationParameter(param)
input  <- "input.demo.vcf"
output <- "out.vcf.gz"
annotateVcf (input, output, param)
pause <- readline("Hit ENTER to continue...")

## ---------- Step 2 ----------
## Now we are able to query the VCF by region or by gene
## Extract by region
##
## A quick way to obtain gentoype matrix:
genotype <- readVCFToMatrixByRange(fileName = output,
                                   range = "1:1-10",
                                   annoType = "")
print(genotype)
pause <- readline("Hit ENTER to continue...")

## You can also fine-control what fields from VCFs to extract:
genotypeList <- readVCFToListByRange(fileName = output,
                                     range = "1:1-10",
                                     annoType = "",
                                     vcfColumn = c("CHROM", "POS"),
                                     vcfInfo = c("ANNO", "SIFT"),
                                     vcfIndv = "GT")
print(genotypeList)
pause <- readline("Hit ENTER to continue...")

## Extract by gene
## Another useful feature for seqminer is to extract data by gene.
## Similar to 'readVCFToMatrixByRange', we use 'readVCFToMatrixByGene':
geneFile = system.file("tabanno/test.gene.txt", package = "seqminer")
genotype <- readVCFToMatrixByGene(fileName = output,
                                  geneFile = geneFile,
                                  geneName = "GENE1",
                                  annoType = "")
print(genotype)
pause <- readline("Hit ENTER to continue...")

## Similar to 'readVCFToMatrixByRange', we use 'readVCFToMatrixByGene':
genotypeList <- readVCFToListByGene(fileName = output,
                                    geneFile = geneFile,
                                    geneName = "GENE1",
                                    annoType = "",
                                    vcfColumn = c("CHROM", "POS"),
                                    vcfInfo = "ANNO",
                                    vcfIndv = "GT")
print(genotypeList)
print(genotypeList)
pause <- readline("Hit ENTER to quit...")

## Thanks for using seqminer
## You can send feedbacks/comments to:
##   Xiaowei Zhan <zhanxw@gmail.com> and
##   Dajiang Liu <dajiang.liu@psu.edu>
## More documentations
##   http://zhanxw.com/seqminer/
