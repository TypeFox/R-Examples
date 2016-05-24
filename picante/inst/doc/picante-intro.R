### R code from vignette source 'picante-intro.Rnw'

###################################################
### code chunk number 1: picante-intro.Rnw:16-17
###################################################
options(width=72)


###################################################
### code chunk number 2: picante-intro.Rnw:30-31
###################################################
library(picante)


###################################################
### code chunk number 3: picante-intro.Rnw:38-43
###################################################
data(phylocom)
names(phylocom)
phy <- phylocom$phylo
comm <- phylocom$sample
traits <- phylocom$traits


###################################################
### code chunk number 4: picante-intro.Rnw:50-51
###################################################
phy


###################################################
### code chunk number 5: picante-intro.Rnw:60-64
###################################################
comm
class(comm)
colnames(comm)
rownames(comm)


###################################################
### code chunk number 6: picante-intro.Rnw:75-78
###################################################
head(traits)
traitA <- df2vec(traits, "traitA")
traitA


###################################################
### code chunk number 7: picante-intro.Rnw:89-91
###################################################
prunedphy <- prune.sample(comm, phy)
prunedphy


###################################################
### code chunk number 8: picante-intro.Rnw:96-101
###################################################
par(mfrow = c(2, 3))
for (i in row.names(comm)) {
plot(prunedphy, show.tip.label = FALSE, main = i)
tiplabels(tip = which(prunedphy$tip.label %in% names(which(comm [i, ] > 0))) , pch=19, cex=2)
}


###################################################
### code chunk number 9: picante-intro.Rnw:106-111
###################################################
par(mfrow=c(2,2))
for (i in names(traits)) {
	plot(phy, show.tip.label=FALSE, main=i)
	tiplabels(pch=22, col=traits[,i]+1, bg=traits[,i]+1, cex=1.5)
}


###################################################
### code chunk number 10: picante-intro.Rnw:122-124
###################################################
pd.result <- pd(comm, phy, include.root=TRUE)
pd.result


###################################################
### code chunk number 11: picante-intro.Rnw:143-148
###################################################
phydist <- cophenetic(phy)
ses.mpd.result <- ses.mpd(comm, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=99)
ses.mpd.result
ses.mntd.result <- ses.mntd(comm, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=99)
ses.mntd.result


###################################################
### code chunk number 12: picante-intro.Rnw:174-175
###################################################
par(mfrow=c(1,1))


###################################################
### code chunk number 13: picante-intro.Rnw:177-181
###################################################
comdist.result <- comdist(comm, phydist)
comdist.result
comdist.clusters <- hclust(comdist.result)
plot(comdist.clusters)


###################################################
### code chunk number 14: picante-intro.Rnw:194-196
###################################################
traits <- traits[phy$tip.label,]
multiPhylosignal(traits,phy)


