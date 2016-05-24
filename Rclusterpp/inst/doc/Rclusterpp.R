### R code from vignette source 'Rclusterpp.Rnw'

###################################################
### code chunk number 1: version
###################################################
prettyVersion <- packageDescription("RcppEigen")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: preliminaries
###################################################
require( Rclusterpp )
Rclusterpp.setThreads(1)


###################################################
### code chunk number 3: simple
###################################################
h <- hclust(dist(USArrests, method="euclidean"), method="average")
r <- Rclusterpp.hclust(USArrests, method="average", distance="euclidean")
# Check equality of the dedrogram tree and agglomeration heights
identical(h$merge, r$merge) && all.equal(h$height, r$height)


###################################################
### code chunk number 4: linkages
###################################################
Rclusterpp.linkageKinds()


###################################################
### code chunk number 5: distances
###################################################
Rclusterpp.distanceKinds()


###################################################
### code chunk number 6: example
###################################################
cat(readLines(system.file("examples","clustering.R",package="Rclusterpp")),sep="\n")


