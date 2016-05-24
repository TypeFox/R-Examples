## This demo will give a typical workflow to read summary statistics
##
##
library(seqminer)

## Now we are able to query the summary statistics file by region or by gene
## e.g. Extract by region
##
## A quick way to obtain all summary statistics from score statistic file.
## Here we used a summary statistics file that have been annotated.
## For how to annotate the summary statistics file, see demo(annotate, package = "seqminer")
score.file <- system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
stats <- rvmeta.readScoreByRange(score.file,
                                    "1:196621000-196623000")
print(names(stats))
print(stats[1:5])

## You can also read covariances from covariance file:
cov.file <- system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
stats <- rvmeta.readCovByRange(cov.file,
                               "1:196621000-196623000")
print(stats)
pause <- readline("Hit ENTER to continue...")

## More commonly, you can read BOTH score statistics and covariance statistics
## e.g. by range
stats <- rvmeta.readDataByRange(score.file, cov.file,
                                "1:196621000-196623000")
print(names(stats))
print(names(stats[[1]]))
print(stats[[1]][1:5])
pause <- readline("Hit ENTER to continue...")
## e.g. by gene
gene.file <- system.file("rvtests/cfh.refFlat.txt", package = "seqminer")
stats <- rvmeta.readDataByGene(score.file, cov.file, geneFile = gene.file,
                               gene = "CFH")

print(names(stats))
print(names(stats[[1]]))
print(stats[[1]][1:5])
pause <- readline("Hit ENTER to quit...")

## Thanks for using seqminer
## You can send feedbacks/comments to:
##   Xiaowei Zhan <zhanxw@gmail.com> and
##   Dajiang Liu <dajiang.liu@psu.edu>
## More documentations
##   http://zhanxw.com/seqminer/
