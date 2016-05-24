### R code from vignette source 'GMD-data-processing.Rnw'

###################################################
### code chunk number 1: GMD-data-processing.Rnw:21-22
###################################################
  cat(as.character(packageVersion('GMD')))


###################################################
### code chunk number 2: GMD-data-processing.Rnw:26-27
###################################################
  cat(unlist(strsplit(packageDescription('GMD')[['Date']],' '))[1])


###################################################
### code chunk number 3: GMD-data-processing.Rnw:196-211
###################################################
require(GMD)

## The file path of ecternal data
id <- "H3K4me1"
inFpath <- system.file('extdata/hg19',paste0(id,'.BedGraph.gz'),package="GMD",mustWork=TRUE)
print(inFpath)

## Convert to a vector of depth-like signals
res <- bedgraph.to.depth(inFpath,chr="chr11",start=17093938,end=17101220,reverse=TRUE)
str(res)

## Visualize the pattern of the signals
plot(as.numeric(names(res)),res,type="l",xlim=rev(range(as.numeric(names(res)))),
     main="H3K4me1 over RPS13", ylab="Depth", xlab="Genomic Postion (chr11)"
     )


###################################################
### code chunk number 4: GMD-data-processing.Rnw:218-226
###################################################
data_test <- list()
ids <- c("H3K4me1","H3K4me2","H3K4me3","H3K9me3","H3K27me3","H3K36me3")
for ( id in ids) {
  inFpath <- system.file('extdata/hg19',paste0(id,'.BedGraph.gz'),package="GMD",mustWork=TRUE) 
  res <- bedgraph.to.depth(inFpath,chr="chr11",start=17093938,end=17101220,reverse=TRUE)
  data_test[[id]] <- res
}
str(data_test)


###################################################
### code chunk number 5: GMD-data-processing.Rnw:233-235
###################################################
x <- gmdm(data_test,sliding=FALSE)
print(x)


###################################################
### code chunk number 6: GMD-data-processing.Rnw:242-243
###################################################
heatmap.3(x,dendrogram="both",cexCol=0.85,cexRow=0.85)


###################################################
### code chunk number 7: GMD-data-processing.Rnw:251-253
###################################################
plot(x,cex.text=2,type="polygon",if.plot.new=FALSE) 
#- setting if.plot.new to TRUE to pop-up a auto-adjust window


