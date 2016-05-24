### R code from vignette source 'vignette1.Rnw'

###################################################
### code chunk number 1: browser1
###################################################
require("rphast")

# extract alignment and annotation files from RPHAST package
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, c("sol1.maf", "sol1.gp"))

# read alignment 
align <- read.msa("sol1.maf")

# read gene annotations from a UCSC "genepred" file
feats <- read.feat("sol1.gp")

# define tree using Newick string
tomatoTree <- "((((tomato, potato), eggplant), pepper), petunia);"


###################################################
### code chunk number 2: browser2
###################################################
names(align)
align$names <- c("tomato", "potato", "pepper", "petunia", "eggplant")
names(align)

unique(feats$seqname)
feats$seqname <- "tomato"


###################################################
### code chunk number 3: browser3
###################################################
feats <- add.introns.feat(feats)
feats <- feats[feats$feature != "exon",]

table(feats$feature)

# make a feature that represents the entire chromosome.  We will ignore 
# several thousand bases at the beginning of the reference genome for 
# which no alignments are available by setting the "start" of the feature 
# equal to the beginning of the aligned region 
wholeChrom <- feat(seq="tomato", src=".", feature="all",
                   start=align$offset, 
                   end=align$offset+ncol.msa(align, "tomato"))

# annotate intergenic regions
intergenicFeats <- inverse.feat(feats, region.bounds=wholeChrom)
intergenicFeats$feature <- "intergenic"
feats <- rbind.feat(feats, intergenicFeats)


###################################################
### code chunk number 4: browser4
###################################################
align4d <- get4d.msa(align, feats)
neutralMod <- phyloFit(align4d, tree=tomatoTree, subst.mod="REV")


###################################################
### code chunk number 5: browser5
###################################################
pc <- phastCons(align, neutralMod, expected.length=6, 
                target.coverage=0.125, viterbi=TRUE)
names(pc)

consElements <- pc$most.conserved
# this shows how many bases are predicted to be conserved
coverage.feat(consElements)

# this shows the fraction of bases covered by conserved elements
coverage.feat(consElements)/coverage.feat(wholeChrom)

# the posterior probabilities for every base are here:
names(pc$post.prob.wig)
dim(pc$post.prob.wig)

# and the overall likelihood is here:
pc$likelihood


###################################################
### code chunk number 6: browser6
###################################################
pp <- phyloP(neutralMod, align, method="LRT", mode="CONACC")

# the returned object is a data frame giving statistics for every base 
# in the alignment
names(pp)
dim(pp)


###################################################
### code chunk number 7: browserFigure
###################################################
codingFeats <- feats[feats$feature=="CDS",]
geneTrack <- as.track.feat(codingFeats, "genes", is.gene=TRUE)
consElTrack <- as.track.feat(consElements, "phastCons most conserved", col="red")
phastConsScoreTrack <- as.track.wig(wig=pc$post.prob.wig,
                                    name="phastCons post prob", col="red", ylim=c(0, 1))
phyloPTrack <- as.track.wig(coord=pp$coord, score=pp$score, name="phyloP score", 
                            col="blue", smooth=TRUE, horiz.line=0)

plot.track(list(geneTrack, consElTrack, phastConsScoreTrack, phyloPTrack),
           xlim=c(60000, 68000), cex.labels=1.25, cex.axis=1.25, cex.lab=1.5)



###################################################
### code chunk number 8: elementLengthByType
###################################################
ce <- pc$most.conserved
plot(density.feat(ce), ylim=c(0, 0.018), 
     main="Element Length by Type", xlab="Length", 
     mgp=c(1.5,0.5,0),mar=c(2,2,2,2))

# obtain elements that overlap codingFeats by at least 50 percent
codingConsEle <- overlap.feat(ce, codingFeats, min.percent=0.5)
# obtain elements that overlap by less than 50 percent
noncodingConsEle <- overlap.feat(ce, codingFeats, min.percent=0.5, 
                            overlapping=FALSE)
lines(density.feat(codingConsEle), col="red")
lines(density.feat(noncodingConsEle), col="blue")
legend(c("All", "Coding", "Noncoding"), x="topright", inset=0.05, 
       lty=1, col=c("black", "red", "blue"))


###################################################
### code chunk number 9: coverage-composition
###################################################
par(mfrow=c(2, 2), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, mar=c(5,5,4,2))

# look at fold-enrichment of each annotation type by conserved element
enrich <- enrichment.feat(ce, feats, wholeChrom)
col <- rainbow(nrow(enrich))
barplot(enrich$enrichment, col=col,
        main="Enrichment of\nConserved Elements",
        ylab="Fold Enrichment")
plot.new()
legend(x="center", legend=enrich$type, fill=col, cex=1.5)
# look at the composition of the conserved elements
comp <- composition.feat(ce, feats)
pie(comp$composition, col=rainbow(nrow(comp)), radius=1.0,
    main="Composition of\nConserved Elements", labels=NA)
# compare with background composition
comp <- composition.feat(wholeChrom, feats)
pie(comp$composition, col=rainbow(nrow(comp)), radius=1.0,
    main="Background\nComposition", labels=NA)


###################################################
### code chunk number 10: em
###################################################
pcEM <- phastCons(align, neutralMod, viterbi=TRUE, estimate.transitions=TRUE) 
names(pcEM)
pcEM$transition.rates
pcEM$likelihood


###################################################
### code chunk number 11: emPlot
###################################################
coverage.feat(pcEM$most.conserved)
coverage.feat(pcEM$most.conserved, pc$most.conserved)
coverage.feat(pcEM$most.conserved, pc$most.conserved, or=TRUE)
coverage.feat(pcEM$most.conserved, pc$most.conserved,
             not=c(FALSE, TRUE))
coverage.feat(pcEM$most.conserved, pc$most.conserved,
             not=c(TRUE, FALSE))
          
plot.track(list(as.track.feat(pc$most.conserved, name="No estimation"),
                as.track.feat(pcEM$most.conserved, name="With estimation")))


###################################################
### code chunk number 12: kdensity
###################################################
plot(density.feat(pc$most.conserved),
     main="Distribution of Element Lengths", xlab="Length", xlim=c(0, 1000))
lines(density.feat(pcEM$most.conserved), col="red")
legend(x="topright", inset=0.05, c("without estimation", "with estimation"), 
       lty=1, col=c("black", "red"))


