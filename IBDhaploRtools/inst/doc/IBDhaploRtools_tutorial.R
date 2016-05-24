### R code from vignette source 'IBDhaploRtools_tutorial.Rnw'

###################################################
### code chunk number 1: load.package
###################################################
## to install the package
## install.packages("IBDhaploRtools", repos="http://R-Forge.R-project.org")

## to load the library at the start of R sesssion
library("IBDhaploRtools")

## to bring up this tutorial
## vignette("IBDhaploRtools_tutorial")


###################################################
### code chunk number 2: load.data
###################################################
## to use the tutorial data
data(qibd_phased)
data(ids_phased)

data(ids_unphased)
data(qibd_unphased)

data(trueibd_phased)

## to use you own data, e.g.
## my.qibd.filename <- '~/Documents/qibd_unphased_2011.gold'
## my.ids.filename <- '~/Documents/ids_unphased_2011.gold'
## my.trueibd.filename <- `~/Documents/outfifteen.txt'

## For data description, e.g.
## ?ids_phased


###################################################
### code chunk number 3: ibdhap.make.calls
###################################################

phased.gold <- ibdhap.make.calls( qibd.file  = qibd_phased,
                   ids.file= ids_phased, cutoff = 0.8)

unphased.gold <- ibdhap.make.calls( qibd.file  = qibd_unphased,
                   ids.file= ids_unphased, cutoff = 0.8)

## or if you specified a file location
## rather than supplying pre-loaded data:
## phased.gold <- ibdhap.make.calls( qibd.filename = my.qibd.filename,
##                    ids.filename = my.ids.filename, cutoff = 0.8)


###################################################
### code chunk number 4: ibdhap.make.calls.show
###################################################
phased.gold[1:20,]


###################################################
### code chunk number 5: ibdhap.names
###################################################
ibdhap.names( ids.file = ids_phased )

  ## or, alternatively
  ## ibdhap.names( ids.filename = my.ids.filename)


###################################################
### code chunk number 6: ibdhap.make.true.show
###################################################
## To load tutorial data:
## data( trueibd_phased )

## To load your own dataset,
## my.trueibd.filename <-`~/Documents/fgl2ibdoutput.txt'
## ibdhap.make.true( true.filename = my.trueibd.filename )


###################################################
### code chunk number 7: ibdhap.summary
###################################################
## summary statistics for phased data:
summary.phased <- ibdhap.summary.calls( phased.gold, data.type="h")
summary.phased

## and for unphased data:
summary.unphased <- ibdhap.summary.calls( unphased.gold, data.type="g")
summary.unphased


###################################################
### code chunk number 8: ibdhap.seg.lengths
###################################################
## load the positions data
data(posvec)

seg.lengths.phased <- ibdhap.seg.lengths(phased.gold[,4], position = posvec)
seg.lengths.phased

seg.lengths.unphased <- ibdhap.seg.lengths(unphased.gold[,4], position = posvec)
seg.lengths.unphased


###################################################
### code chunk number 9: ibdhap.transitions
###################################################
##set the display so the matrix displays nicely
options(width=100)

## transition matrix
transitions.phased <- ibdhap.transitions(phased.gold,
                                         data.type="h")
transitions.unphased <- ibdhap.transitions(unphased.gold,
                                           data.type="g")

## print the transitions matrix into two parts so it
## fits on the document
transitions.phased[,1:8]
transitions.phased[,9:15]

## unphased transitions matrix
transitions.unphased


###################################################
### code chunk number 10: ibdhap.barplot.phased
###################################################
## Figure 1
par(mfrow=c(4,1))
ibdhap.barplot(phased.gold[,1], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(phased.gold[,2], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(phased.gold[,3], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(phased.gold[,4], data.type="h",
               xlab="", ylab="", position = posvec)


###################################################
### code chunk number 11: ibdhap.barplot.unphased
###################################################
## Figure 2
par(mfrow=c(4,1))
ibdhap.barplot(unphased.gold[,1], data.type="g",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(unphased.gold[,2], data.type="g",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(unphased.gold[,3], data.type="g",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(unphased.gold[,4], data.type="g",
               xlab="", ylab="", position = posvec)


###################################################
### code chunk number 12: fig1
###################################################
## Figure 1
par(mfrow=c(4,1))
ibdhap.barplot(phased.gold[,1], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(phased.gold[,2], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(phased.gold[,3], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(phased.gold[,4], data.type="h",
               xlab="", ylab="", position = posvec)


###################################################
### code chunk number 13: fig2
###################################################
## Figure 2
par(mfrow=c(4,1))
ibdhap.barplot(unphased.gold[,1], data.type="g",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(unphased.gold[,2], data.type="g",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(unphased.gold[,3], data.type="g",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(unphased.gold[,4], data.type="g",
               xlab="", ylab="", position = posvec)


###################################################
### code chunk number 14: load.simulated.gold
###################################################
## load the true IBD states data set
data(trueibd_phased)

## map phased to unphased labels
trueibd_unphased <- h.to.g( trueibd_phased )

## view the first few rows of each data set
head(trueibd_phased)
head(trueibd_unphased)


###################################################
### code chunk number 15: trueibd.phased.barplot
###################################################
## Figure 3
par(mfrow=c(4,1))
ibdhap.barplot(trueibd_phased[,1], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(trueibd_phased[,2], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(trueibd_phased[,3], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(trueibd_phased[,4], data.type="h",
               xlab="", ylab="", position = posvec)


###################################################
### code chunk number 16: fig3
###################################################
## Figure 3
par(mfrow=c(4,1))
ibdhap.barplot(trueibd_phased[,1], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(trueibd_phased[,2], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(trueibd_phased[,3], data.type="h",
               xlab="", ylab="", position = posvec)
ibdhap.barplot(trueibd_phased[,4], data.type="h",
               xlab="", ylab="", position = posvec)


###################################################
### code chunk number 17: ibd.compare.loci.show
###################################################
## phased data
ibdhap.compare.loci(phased.gold, trueibd_phased, "h")
## unphased data
ibdhap.compare.loci(unphased.gold, trueibd_unphased, "g")


###################################################
### code chunk number 18: ibdhap.compare.segs.show
###################################################
## Run the function on the phased data
## NB: we can also use pos = posvec
corr.seg.phased <- ibdhap.compare.segs( phased.gold,
                                       trueibd_phased,
                                       "h", 0.8, pos = NA)

## First element of output is summary statistics
corr.seg.phased$seg.stats

## Same thing for the unphased data
## first converting true states to unphased
trueibd_unphased <- h.to.g( trueibd_phased )
corr.seg.unphased <- ibdhap.compare.segs( unphased.gold,
                                         trueibd_unphased,
                                         "g", 0.8, pos = NA )
corr.seg.unphased$seg.stats



###################################################
### code chunk number 19: fig4.code
###################################################
par(mfrow=c(1,1)) ## reset the plot window to square
## Figure 4
plot( corr.seg.phased$seg.info[,1], corr.seg.phased$seg.info[,5],
       xlab="segment length", ylab="proportion correct",
       main = "Phased Data")

 lines(lowess(corr.seg.phased$seg.info[,1],
       corr.seg.phased$seg.info[,5]), col="blue")


###################################################
### code chunk number 20: fig5.code
###################################################
## Figure 5
plot( corr.seg.unphased$seg.info[,1], corr.seg.unphased$seg.info[,5],
       xlab="segment length", ylab="proportion correct",
       main = "Unphased Data")

lines(lowess(corr.seg.unphased$seg.info[,1],
      corr.seg.unphased$seg.info[,5]), col="blue")


###################################################
### code chunk number 21: fig4
###################################################
par(mfrow=c(1,1)) ## reset the plot window to square
## Figure 4
plot( corr.seg.phased$seg.info[,1], corr.seg.phased$seg.info[,5],
       xlab="segment length", ylab="proportion correct",
       main = "Phased Data")

 lines(lowess(corr.seg.phased$seg.info[,1],
       corr.seg.phased$seg.info[,5]), col="blue")


###################################################
### code chunk number 22: fig5
###################################################
## Figure 5
plot( corr.seg.unphased$seg.info[,1], corr.seg.unphased$seg.info[,5],
       xlab="segment length", ylab="proportion correct",
       main = "Unphased Data")

lines(lowess(corr.seg.unphased$seg.info[,1],
      corr.seg.unphased$seg.info[,5]), col="blue")


