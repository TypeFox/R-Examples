### R code from vignette source 'SigTree.Rnw'

###################################################
### code chunk number 1: SigTree.Rnw:102-103 (eval = FALSE)
###################################################
## tre.path <- "C:/folder1/folder2"


###################################################
### code chunk number 2: SigTree.Rnw:109-110
###################################################
tre.path <- system.file("sample", package="SigTree")


###################################################
### code chunk number 3: SigTree.Rnw:121-125
###################################################
library(ape)
tree.file <- paste(tre.path,"sample.tre",sep="/")
tree <- read.tree(tree.file)
plot(tree,type="fan")


###################################################
### code chunk number 4: SigTree.Rnw:145-148
###################################################
sig.file <- paste(tre.path,"sample.csv",sep="/")
frame <- read.csv(sig.file)
head(frame)


###################################################
### code chunk number 5: SigTree.Rnw:163-164
###################################################
library(SigTree)


###################################################
### code chunk number 6: SigTree.Rnw:177-178
###################################################
plotSigTree(tree, frame)


###################################################
### code chunk number 7: SigTree.Rnw:215-220
###################################################
library(RColorBrewer)
RdBu <- brewer.pal(7, "RdBu")
RdBu[4] <- brewer.pal(7, "Greys")[3]
plotSigTree(tree, frame, pal=RdBu, tip.label.size=.8, edge.width=2)
legend("topright",c("C > T","C < T"),cex=1.5,text.col=RdBu[c(7,1)],bty="n")


###################################################
### code chunk number 8: SigTree.Rnw:243-244 (eval = FALSE)
###################################################
## plotSigTree(tree, frame, side=2)


###################################################
### code chunk number 9: SigTree.Rnw:253-257
###################################################
plotSigTree(tree, frame, pal=RdBu, tip.label.size=.8, edge.width=2)
p.cut.leg <- c("0.00 - 0.01","0.01 - 0.05","0.05 - 0.10",
               "0.10 - 0.90","0.90 - 0.95","0.95 - 0.99","0.99 - 1.00")
legend("topright",rev(p.cut.leg),text.col=rev(RdBu), bty="n", cex=1)


###################################################
### code chunk number 10: SigTree.Rnw:297-302
###################################################
plotSigTree(tree, frame, pal=RdBu, tip.label.size=.8, edge.width=2, 
            branch="node")
p.cut.leg <- c("0.00 - 0.01","0.01 - 0.05","0.05 - 0.10",
               "0.10 - 0.90","0.90 - 0.95","0.95 - 0.99","0.99 - 1.00")
legend("topright",rev(p.cut.leg),text.col=rev(RdBu), bty="n", cex=1)


###################################################
### code chunk number 11: SigTree.Rnw:318-319
###################################################
p.adjust.methods


###################################################
### code chunk number 12: SigTree.Rnw:328-329
###################################################
plotSigTree(tree, frame, method="BY")


###################################################
### code chunk number 13: SigTree.Rnw:345-348
###################################################
plotSigTree(tree, frame, pal=RdBu, tip.label.size=.8, edge.width=2, 
            branch.label=TRUE, branch.label.size=.75)
edgelabels(edge=176,frame="circ",bg="yellow",cex=.8)


###################################################
### code chunk number 14: SigTree.Rnw:356-357 (eval = FALSE)
###################################################
## export.inherit(tree, frame, file="sampleInherit.csv")


###################################################
### code chunk number 15: SigTree.Rnw:362-363
###################################################
temp <- export.inherit(tree, frame, frame=TRUE)


###################################################
### code chunk number 16: SigTree.Rnw:369-372
###################################################
br176 <- temp[temp$Branch=="176",]
t <- !is.na(br176[1,])
br176[1,t]


###################################################
### code chunk number 17: SigTree.Rnw:380-383
###################################################
br106 <- temp[temp$Branch=="106",]
t <- !is.na(br106[1,])
br106[1,t]


###################################################
### code chunk number 18: SigTree.Rnw:393-395
###################################################
library(phyext2) 
export.figtree(tree, frame, pal=RdBu, file="sigsample.tre")


###################################################
### code chunk number 19: SigTree.Rnw:420-424
###################################################
keep.taxa <- c("t57","t99","t53","t62","t39","t16",
               "t63","t67","t45","t1","t34","t82")
library(phyloseq)
new_tree <- prune_taxa(keep.taxa, tree)


###################################################
### code chunk number 20: SigTree.Rnw:430-434
###################################################
t <- is.element(frame$OTU, keep.taxa)
new_frame <- frame[t,]
plotSigTree(new_tree, new_frame, pal=RdBu, tip.label.size=1.5, 
            edge.width=4)


###################################################
### code chunk number 21: SigTree.Rnw:472-473 (eval = FALSE)
###################################################
## sig_tree <- read.tree("sigsample.tre")


###################################################
### code chunk number 22: SigTree.Rnw:486-487
###################################################
sig_tree <- read.nexus("sigsample.tre")


###################################################
### code chunk number 23: SigTree.Rnw:495-497 (eval = FALSE)
###################################################
## singletontree.file <- paste(tre.path,"singletonsample.tre",sep="/")
## tree <- read.tree(singletontree.file)


###################################################
### code chunk number 24: SigTree.Rnw:515-528 (eval = FALSE)
###################################################
## t1 <- read.table(singletontree.file)
## # remove initial (
## t1.split <- strsplit(as.character(t1$V1),"")[[1]]
## t2 <- paste(t1.split[-1],collapse="")
## # remove last ) and the num:num preceding it
## t2.split <- strsplit(t2,")")[[1]]
## t2.len <- length(t2.split)
## t3.split <- t2.split[c(1:(t2.len-2),t2.len)]
## t3 <- paste(t3.split,collapse=")")
## # write to file (fixed now) and read in again
## newtree.file <- paste(tre.path,"fixedsinglesample.tre",sep="/")
## write.table(t3,file=newtree.file, quote=F, col.names=F, row.names=F)
## t4 <- read.tree(newtree.file)


