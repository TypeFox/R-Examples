### R code from vignette source 'PlotsAndStats.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PlotsAndStats.Rnw:33-46
###################################################
options(warn=2)

library(cheddar)

# Makes copy-paste much less painful
options(continue=' ')
options(width=90)
options(prompt='> ')

options(SweaveHooks = list(fig=function() par(mgp=c(2.5,1,0), 
                                              mar=c(4,4,2,1),
                                              oma=c(0,0,1,0),
                                              cex.main=0.8)))


###################################################
### code chunk number 2: PlotsAndStats.Rnw:124-126
###################################################
getOption("SweaveHooks")[["fig"]]()
data(TL84)
PlotNPS(TL84, 'Log10M', 'Log10N')


###################################################
### code chunk number 3: PlotsAndStats.Rnw:144-145
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPS(TL84, 'Log10M', 'Log10N', show.web=FALSE, highlight.nodes=NULL)


###################################################
### code chunk number 4: PlotsAndStats.Rnw:152-153
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPS(TL84, 'Log10M', 'Log10N', show.nodes.as='labels', show.web=FALSE)


###################################################
### code chunk number 5: PlotsAndStats.Rnw:160-162
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPS(TL84, 'Log10M', 'Log10N', show.nodes.as='labels', show.web=FALSE, 
        node.labels='node', cex=0.5)


###################################################
### code chunk number 6: PlotsAndStats.Rnw:168-171
###################################################
getOption("SweaveHooks")[["fig"]]()
lots.of.letters <- c(letters, LETTERS, paste(LETTERS,letters,sep=''))
PlotNPS(TL84, 'Log10M', 'Log10N', show.nodes.as='labels', show.web=FALSE, 
        node.labels=lots.of.letters[1:NumberOfNodes(TL84)])


###################################################
### code chunk number 7: PlotsAndStats.Rnw:176-177
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPS(TL84, 'Log10M', 'Log10N', show.nodes.as='both', show.web=FALSE, cex=2)


###################################################
### code chunk number 8: PlotsAndStats.Rnw:187-189
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPS(TL84, 'Log10M', 'Log10N', xlab=Log10MLabel(TL84), 
        ylab=Log10NLabel(TL84))


###################################################
### code chunk number 9: PlotsAndStats.Rnw:199-208
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,3))
PlotNPS(TL84, 'Log10M', 'OutDegree', show.web=FALSE)
abline(lm(OutDegree(TL84) ~ Log10M(TL84)))

PlotNPS(TL84, 'Log10M', 'InDegree', show.web=FALSE)
abline(lm(InDegree(TL84) ~ Log10M(TL84)))

PlotNPS(TL84, 'Log10M', 'Degree', show.web=FALSE)
abline(lm(Degree(TL84) ~ Log10M(TL84)))


###################################################
### code chunk number 10: PlotsAndStats.Rnw:220-221
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPS(TL84, 'Log10M', 'PreyAveragedTrophicLevel')


###################################################
### code chunk number 11: PlotsAndStats.Rnw:228-229
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPS(TL84, 'Log10M', 'ChainAveragedTrophicLevel')


###################################################
### code chunk number 12: PlotsAndStats.Rnw:244-249
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,2))
PlotNPS(TL84, 'Log10M', 'PreyAveragedTrophicLevel', ylim=c(1, 6), 
        main='Prey-averaged')
PlotNPS(TL84, 'Log10M', 'ChainAveragedTrophicLevel', ylim=c(1, 6), 
        main='Chain-averaged')


###################################################
### code chunk number 13: PlotsAndStats.Rnw:263-268
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
PlotMvN(TL84)
PlotNvM(TL84)
PlotBvM(TL84)
PlotMvB(TL84)


###################################################
### code chunk number 14: PlotsAndStats.Rnw:281-282
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotRankNPS(TL84, 'Log10N')


###################################################
### code chunk number 15: PlotsAndStats.Rnw:287-288
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotRankNPS(TL84, 'Log10N', rank.by='M')


###################################################
### code chunk number 16: PlotsAndStats.Rnw:296-297
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotRankNPS(TL84, 'Log10N', rank.by='M', show.web=TRUE)


###################################################
### code chunk number 17: PlotsAndStats.Rnw:302-303
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotRankNPS(TL84, 'PreyAveragedTrophicLevel', rank.by='M')


###################################################
### code chunk number 18: PlotsAndStats.Rnw:311-312
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotRankNPS(TL84, 'PreyAveragedTrophicLevel', rank.by='M', log10.rank=TRUE)


###################################################
### code chunk number 19: PlotsAndStats.Rnw:322-326
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,3))
PlotMvRankM(TL84)
PlotNvRankN(TL84)
PlotBvRankB(TL84)


###################################################
### code chunk number 20: PlotsAndStats.Rnw:340-341
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPSDistribution(TL84, 'Log10M')


###################################################
### code chunk number 21: PlotsAndStats.Rnw:347-348
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNPSDistribution(TL84, 'Log10M', density.args=list(bw=3))


###################################################
### code chunk number 22: PlotsAndStats.Rnw:368-369
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, col=1, pch=19, highlight.nodes=NULL)


###################################################
### code chunk number 23: PlotsAndStats.Rnw:376-377
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, col=1:56, pch=19, highlight.nodes=NULL)


###################################################
### code chunk number 24: PlotsAndStats.Rnw:386-387
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, colour.by='resolved.to', pch=19, highlight.nodes=NULL)


###################################################
### code chunk number 25: PlotsAndStats.Rnw:395-399
###################################################
getOption("SweaveHooks")[["fig"]]()
colour.spec <- c(Species='purple3', Genus='green3', 'red3')
PlotNvM(TL84, colour.by='resolved.to', colour.spec=colour.spec, pch=19, 
        highlight.nodes=NULL)
legend("topright", legend=names(colour.spec), pch=19, col=colour.spec)


###################################################
### code chunk number 26: PlotsAndStats.Rnw:410-422
###################################################
getOption("SweaveHooks")[["fig"]]()
symbol.spec = c(Bacteria=21, Plantae=22, Chromista=23, 
                Protozoa=24, Animalia=25, 19)
colour.spec = c(Bacteria='purple3', Plantae='green3', 
                Chromista='blue3', Protozoa='orange3', 
                Animalia='red3', 'black')
PlotNvM(TL84, 
        symbol.by='kingdom', symbol.spec=symbol.spec, 
        bg.by='kingdom', bg.spec=colour.spec, 
        colour.by='kingdom', colour.spec=colour.spec, 
        highlight.nodes=NULL)
legend("topright", legend=names(colour.spec), pch=symbol.spec, 
       col=colour.spec, pt.bg=colour.spec)


###################################################
### code chunk number 27: PlotsAndStats.Rnw:434-448
###################################################
getOption("SweaveHooks")[["fig"]]()
symbol.spec = c(Bacteria=21, Plantae=22, Chromista=23, 
                Protozoa=24, Animalia=25, 19)
colour.spec = c(Bacteria='purple3', Plantae='green3', 
                Chromista='blue3', Protozoa='orange3', 
                Animalia='red3', 'black')
PlotNvM(TL84, 
        symbol.by='kingdom', symbol.spec=symbol.spec, 
        bg.by='kingdom', bg.spec=colour.spec, 
        colour.by='kingdom', colour.spec=colour.spec, 
        highlight.nodes=NULL, show.web=FALSE)
legend("topright", legend=names(colour.spec), pch=symbol.spec, 
       col=colour.spec, pt.bg=colour.spec)
models <- NvMLinearRegressions(TL84, class='kingdom')
colours <- PlotLinearModels(models, colour.spec=colour.spec)


###################################################
### code chunk number 28: PlotsAndStats.Rnw:459-460
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, pch=NA, highlight.nodes=NULL)


###################################################
### code chunk number 29: PlotsAndStats.Rnw:473-482
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,2))

# Don't add ticks
options(cheddarTopAndRightTicks=FALSE)
PlotNvM(TL84)

# Add ticks
options(cheddarTopAndRightTicks=TRUE)
PlotNvM(TL84)


###################################################
### code chunk number 30: PlotsAndStats.Rnw:498-499
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, highlight.nodes=Cannibals)


###################################################
### code chunk number 31: PlotsAndStats.Rnw:505-506
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, highlight.nodes=IsolatedNodes)


###################################################
### code chunk number 32: PlotsAndStats.Rnw:512-513
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, highlight.nodes='Chaoborus punctipennis')


###################################################
### code chunk number 33: PlotsAndStats.Rnw:526-527
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, highlight.links=ResourceLargerThanConsumer)


###################################################
### code chunk number 34: PlotsAndStats.Rnw:533-535
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(TL84, highlight.nodes='Chaoborus punctipennis', 
        highlight.links=TrophicLinksForNodes(TL84, 'Chaoborus punctipennis'))


###################################################
### code chunk number 35: PlotsAndStats.Rnw:556-560
###################################################
getOption("SweaveHooks")[["fig"]]()
data(YthanEstuary)
par(mfrow=c(1,2))
PlotNvM(YthanEstuary)
PlotNvM(YthanEstuary, show.na=TRUE)


###################################################
### code chunk number 36: PlotsAndStats.Rnw:570-571
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotNvM(YthanEstuary, xlim=c(-10, 4), ylim=c(-10, 13), show.na=TRUE)


###################################################
### code chunk number 37: PlotsAndStats.Rnw:583-606
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
np <- NPS(TL84)
np[1,'M'] <- NA
PlotNvM(Community(nodes=np, trophic.links=TLPS(TL84), properties=CPS(TL84)), 
        main='Node 1 M=NA', show.nodes.as='both', cex=2, show.na=TRUE)

np <- NPS(TL84)
np[1,'N'] <- NA
PlotNvM(Community(nodes=np, trophic.links=TLPS(TL84), properties=CPS(TL84)), 
        main='Node 1 N=NA', show.nodes.as='both', cex=2, show.na=TRUE)

np <- NPS(TL84)
np[1,'M'] <- NA
np[1,'N'] <- NA
PlotNvM(Community(nodes=np, trophic.links=TLPS(TL84), properties=CPS(TL84)), 
        main='Node 1 M=NA and N=NA', show.nodes.as='both', cex=2, show.na=TRUE)

np <- NPS(TL84)
np[c(10, 20, 30, 40),'M'] <- NA
np[c(10, 20, 30, 40),'N'] <- NA
PlotNvM(Community(nodes=np, trophic.links=TLPS(TL84), properties=CPS(TL84)), 
        main='Nodes 10, 20, 30 and 40 M=NA and N=NA', show.nodes.as='both', 
        cex=2, show.na=TRUE)


###################################################
### code chunk number 38: PlotsAndStats.Rnw:618-622
###################################################
getOption("SweaveHooks")[["fig"]]()
data(YthanEstuary)
par(mfrow=c(1,2))
PlotMvRankM(YthanEstuary)
PlotMvRankM(YthanEstuary, show.na=TRUE)


###################################################
### code chunk number 39: PlotsAndStats.Rnw:646-647
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotTLPS(TL84, 'resource.Log10M', 'consumer.Log10M')


###################################################
### code chunk number 40: PlotsAndStats.Rnw:656-657
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotTLPS(TL84, 'resource.Log10M', 'consumer.Log10M', axes.limits.equal=TRUE)


###################################################
### code chunk number 41: PlotsAndStats.Rnw:678-683
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
PlotPredationMatrix(TL84)
PlotMRvMC(TL84)
PlotNCvNR(TL84)
PlotBRvBC(TL84)


###################################################
### code chunk number 42: PlotsAndStats.Rnw:697-698
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotMRvMC(TL84)


###################################################
### code chunk number 43: PlotsAndStats.Rnw:706-708
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotMRvMC(TL84, colour.by='consumer.category', bg.by='consumer.category', 
          symbol.by='consumer.category')


###################################################
### code chunk number 44: PlotsAndStats.Rnw:719-722
###################################################
SumMByClass(TL84)
SumNByClass(TL84)
SumBiomassByClass(TL84)


###################################################
### code chunk number 45: PlotsAndStats.Rnw:726-729
###################################################
SumMByClass(TL84, 'kingdom')
SumNByClass(TL84, 'kingdom')
SumBiomassByClass(TL84, 'kingdom')


###################################################
### code chunk number 46: PlotsAndStats.Rnw:735-737
###################################################
SumBiomassByClass(TL84)
ApplyByClass(TL84, 'Biomass', 'category', sum)


###################################################
### code chunk number 47: PlotsAndStats.Rnw:750-752
###################################################
models <- NvMLinearRegressions(TL84)
names(models)


###################################################
### code chunk number 48: PlotsAndStats.Rnw:755-756
###################################################
sapply(models, 'coef')


###################################################
### code chunk number 49: PlotsAndStats.Rnw:763-765
###################################################
models <- NvMLinearRegressions(TL84, class='phylum')
names(models)


###################################################
### code chunk number 50: PlotsAndStats.Rnw:773-774
###################################################
sapply(models, is.null)


###################################################
### code chunk number 51: PlotsAndStats.Rnw:780-783
###################################################
data(BroadstoneStream)
models <- NvMLinearRegressions(BroadstoneStream)
sapply(models, is.null)


###################################################
### code chunk number 52: PlotsAndStats.Rnw:787-790
###################################################
NvMSlope(TL84)
NvMIntercept(TL84)
NvMSlopeAndIntercept(TL84)


###################################################
### code chunk number 53: PlotsAndStats.Rnw:793-796
###################################################
NvMSlopeByClass(TL84)
NvMInterceptByClass(TL84)
NvMSlopeAndInterceptByClass(TL84)


###################################################
### code chunk number 54: PlotsAndStats.Rnw:799-802
###################################################
NvMSlopeByClass(TL84, class='kingdom')
NvMInterceptByClass(TL84, class='kingdom')
NvMSlopeAndInterceptByClass(TL84, class='kingdom')


###################################################
### code chunk number 55: PlotsAndStats.Rnw:838-845
###################################################
getOption("SweaveHooks")[["fig"]]()
data(TL84, TL86)
par(mfrow=c(1,2))
PlotMvN(TL84, show.nodes.as='both', cex=2, xlim=c(-2, 10), ylim=c(-14, 0), 
        highlight.nodes=NULL, highlight.links=NULL, main='')
PlotMvN(TL86, show.nodes.as='both', cex=2, xlim=c(-2, 10), ylim=c(-14, 0), 
        highlight.nodes=NULL, highlight.links=NULL, main='')
title(main='Jonsson et al. (2005) AER, Fig. 3 (p 30)', outer=TRUE)


###################################################
### code chunk number 56: PlotsAndStats.Rnw:854-860
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,2))
PlotMCvMR(TL84, xlim=c(-14, 0), ylim=c(-14, 0), main='')
abline(a=0, b=1, lty=2)
PlotMCvMR(TL86, xlim=c(-14, 0), ylim=c(-14, 0), main='')
abline(a=0, b=1, lty=2)
title(main='Jonsson et al. (2005) AER, Fig. 4 (p 33)', outer=TRUE)


###################################################
### code chunk number 57: PlotsAndStats.Rnw:869-875
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
PlotNvM(TL84, xlim=c(-14, 0), ylim=c(-2,10), show.web=FALSE, main='')
PlotNvM(TL86, xlim=c(-14, 0), ylim=c(-2,10), show.web=FALSE, main='')
PlotBvM(TL84, xlim=c(-14, 0), ylim=c(-8,2), show.web=FALSE, main='')
PlotBvM(TL86, xlim=c(-14, 0), ylim=c(-8,2), show.web=FALSE, main='')
title(main='Jonsson et al. (2005) AER, Fig. 5 (p 37)', outer=TRUE)


###################################################
### code chunk number 58: PlotsAndStats.Rnw:884-894
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
PlotNCvNR(TL84, xlim=c(0, 10), ylim=c(-2,10), main='')
abline(a=0, b=1, lty=2)
PlotNCvNR(TL86, xlim=c(0, 10), ylim=c(-2,10), main='')
abline(a=0, b=1, lty=2)
PlotBCvBR(TL84, xlim=c(-8, -2), ylim=c(-8, -2), main='')
abline(a=0, b=1, lty=2)
PlotBCvBR(TL86, xlim=c(-8, -2), ylim=c(-8, -2), main='')
abline(a=0, b=1, lty=2)
title(main='Jonsson et al. (2005) AER, Fig. 7 (p 47)', outer=TRUE)


###################################################
### code chunk number 59: PlotsAndStats.Rnw:903-913
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
TL84.no.iso <- RemoveIsolatedNodes(TL84)
TL86.no.iso <- RemoveIsolatedNodes(TL86)
tl84.levels <- floor(TrophicHeight(TL84.no.iso))
tl86.levels <- floor(TrophicHeight(TL86.no.iso))
PlotNPyramid(TL84.no.iso, level=tl84.levels, main='', ylab='Trophic height')
PlotNPyramid(TL86.no.iso, level=tl86.levels, main='')
PlotBPyramid(TL84.no.iso, level=tl84.levels, main='', ylab='Trophic height')
PlotBPyramid(TL86.no.iso, level=tl86.levels, main='')
title(main='Jonsson et al. (2005) AER, Fig. 8 (p 49)', outer=TRUE)


###################################################
### code chunk number 60: PlotsAndStats.Rnw:922-928
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
PlotNvRankN(TL84, xlim=c(0,60), ylim=c(-2, 10), main='')
PlotNvRankN(TL86, xlim=c(0,60), ylim=c(-2, 10), main='')
PlotBvRankB(TL84, xlim=c(0,60), ylim=c(-8, -2), main='')
PlotBvRankB(TL86, xlim=c(0,60), ylim=c(-8, -2), main='')
title(main='Jonsson et al. (2005) AER, Fig. 10 (p 57)', outer=TRUE)


###################################################
### code chunk number 61: PlotsAndStats.Rnw:937-949
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
PlotRankNPS(TL84, property='Log10N', rank.by='M', log10.rank=TRUE, 
            xlim=c(0,2), ylim=c(-2, 10), ylab=Log10NLabel(TL84), main='')
PlotRankNPS(TL86, property='Log10N', rank.by='M', log10.rank=TRUE, 
            xlim=c(0,2), ylim=c(-2, 10), ylab=Log10NLabel(TL84), main='')
PlotRankNPS(TL84, property='Log10Biomass', rank.by='M', 
            log10.rank=TRUE, xlim=c(0,2), ylim=c(-8, -2), 
            ylab=Log10BLabel(TL84), main='')
PlotRankNPS(TL86, property='Log10Biomass', rank.by='M', 
            log10.rank=TRUE, xlim=c(0,2), ylim=c(-8, -2), 
            ylab=Log10BLabel(TL84), main='')
title(main='Jonsson et al. (2005) AER, Fig. 11 (p 60)', outer=TRUE)


###################################################
### code chunk number 62: PlotsAndStats.Rnw:960-984
###################################################
getOption("SweaveHooks")[["fig"]]()
PlotCommunityVCommunity <- function(a, b, property, xlim=NULL, ylim=NULL, ...)
{
    a.nodes <- NP(a, 'node')
    b.nodes <- NP(b, 'node')
    all.nodes <- union(a.nodes, b.nodes)

    a.values <- NPS(a, property)[,property]
    names(a.values) <- a.nodes
    b.values <- NPS(b, property)[,property]
    names(b.values) <- b.nodes
    points <- PlaceMissingPoints(a.values[all.nodes], xlim,
                                 b.values[all.nodes], ylim)
    plot(points[,1], points[,2], xlim=xlim, ylim=ylim, ...)

    abline(a=0, b=1, lty=2)
}

par(mfrow=c(1,2))
PlotCommunityVCommunity(TL84, TL86, 'Log10N', xlim=c(-2,10), ylim=c(-2,10), 
                        xlab=~log[10]~(N~of~84), ylab=~log[10]~(N~of~86),pch=19)
PlotCommunityVCommunity(TL84, TL86, 'Log10Biomass', 
                        xlim=c(-8,-2), ylim=c(-8,-2), 
                        xlab=~log[10]~(B~of~84), ylab=~log[10]~(B~of~86),pch=19)
title(main='Jonsson et al. (2005) AER, Fig. 12 (p 61)', outer=TRUE)


###################################################
### code chunk number 63: PlotsAndStats.Rnw:997-1013
###################################################
getOption("SweaveHooks")[["fig"]]()
data(pHWebs)
par(mfrow=c(2,2))
for(community in pHWebs[1:2])
{
    PlotNvM(community, xlim=c(-15, 10), ylim=c(-5,15), main='', 
            highlight.nodes=NULL)
    text(-15, 13, with(CPS(community), paste(title, ', pH ', pH, sep='')), 
         adj=0, cex=1.5)
    tlps <- TLPS(community, node.properties='M')
    tlps <- tlps[!is.na(tlps$resource.M) & !is.na(tlps$consumer.M),]
    interaction.strength <- log10( (tlps$consumer.M / tlps$resource.M)^0.75 )
    plot(density(interaction.strength), xlim=c(-4,14), ylim=c(0,0.6), 
         main='', xlab=~log[10]((M[C]/M[R])^0.75))
    rug(interaction.strength)
}
title(main='Layer et al. (2010) AER, Fig. 6 (p 282)', outer=TRUE)


###################################################
### code chunk number 64: PlotsAndStats.Rnw:1025-1040
###################################################
getOption("SweaveHooks")[["fig"]]()
data(BroadstoneStream)
par(mfrow=c(1,2))
PlotMvN(BroadstoneStream, show.nodes.as='labels', label.cex=0.8, 
        xlim=c(-2, 4.2), ylim=c(-6,2), main='', show.na=FALSE, 
        highlight.links=NULL)
abline(a=0, b=-1)

tlps <- TLPS(BroadstoneStream, node.properties='M')
lty <- rep(0, NumberOfTrophicLinks(BroadstoneStream))
lty[tlps$resource.M > tlps$consumer.M] <- 1
PlotMvN(BroadstoneStream, show.nodes.as='labels', label.cex=0.8, 
        xlim=c(-2, 4.2), ylim=c(-6,2), main='', show.na=FALSE, 
        highlight.links=NULL, link.lty=lty)
abline(a=0, b=-1)
title(main='Woodward et al. (2005) AER, Fig. 4 (p 108)', outer=TRUE)


###################################################
### code chunk number 65: PlotsAndStats.Rnw:1055-1059
###################################################
data(TL84, TL86, YthanEstuary)
collection <- CommunityCollection(list(TL84, TL86, YthanEstuary))
table <- NvMTriTrophicTable(collection)
print(round(table,2))


###################################################
### code chunk number 66: PlotsAndStats.Rnw:1063-1092
###################################################
res <- lapply(list(TL84, TL86, YthanEstuary), function(community)
{
    community <- RemoveNodes(community, remove=with(NPS(community), node[is.na(M) | is.na(N)]))
    community <- RemoveCannibalisticLinks(community)
    community <- RemoveIsolatedNodes(community)

    chains <- ThreeNodeChains(community, node.properties='M')
    MR <- chains$bottom.M
    MI <- chains$intermediate.M
    MC <- chains$top.M

    lp <- TLPS(community, node.properties='M')

    return (c('MR<=MI<=MC'=sum(MR<=MI & MI<=MC), 
              'MR<=MC<MI'=sum(MR<=MC & MC<MI), 
              'MI<MR<=MC'=sum(MI<MR  & MR<=MC), 
              'MI<=MC<MR'=sum(MI<=MC & MC<MR), 
              'MC<MR<MI'=sum(MC<MR  & MR<MI), 
              'MC<MI<MR'=sum(MC<MI  & MI<MR), 
              'All 2-chains'=nrow(chains),
              'MR<MC'=sum(lp$resource.M<lp$consumer.M), 
              'MR=MC'=sum(lp$resource.M==lp$consumer.M), 
              'MR>MC'=sum(lp$resource.M>lp$consumer.M), 
              'All links'=nrow(lp)))
})
res <- do.call('cbind', res)

colnames(res) <- c('TL84', 'TL86', 'Ythan Estuary')
print(round(res,2))


###################################################
### code chunk number 67: PlotsAndStats.Rnw:1099-1112
###################################################
getOption("SweaveHooks")[["fig"]]()
data(TL84, TL86, YthanEstuary)
par(mfrow=c(3,2))
for(community in list(TL84, TL86, YthanEstuary))
{
    community <- RemoveIsolatedNodes(community)
    pch <- rep(1, NumberOfNodes(community))
    pch[IsIntermediateNode(community)] <- 20
    pch[IsTopLevelNode(community)] <- 8
    PlotNvM(community, col=1, highlight.nodes=NULL, show.web=FALSE, 
            main='', pch=pch)
    PlotAuppervAlower(community, main='')
}
title(main='Cohen et al. (2009) PNAS, Fig. 1 (p 22336)', outer=TRUE)


###################################################
### code chunk number 68: PlotsAndStats.Rnw:1121-1124
###################################################
data(ChesapeakeBay)
res <- NodeQuantitativeDescriptors(ChesapeakeBay, 'biomass.flow')
print(round(res[1:6,],2))


###################################################
### code chunk number 69: PlotsAndStats.Rnw:1128-1131
###################################################
data(ChesapeakeBay)
res <- QuantitativeDescriptors(ChesapeakeBay, 'biomass.flow')
print(round(res,3))


