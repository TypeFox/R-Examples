### R code from vignette source 'Collections.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Collections.Rnw:32-45
###################################################
options(warn=2)

library(cheddar)

# Makes copy-paste much less painful
options(continue=' ')
options(prompt='> ')
options(width=90)

options(SweaveHooks = list(fig=function() par(mgp=c(2.5,1,0), 
                                              mar=c(4,4,2,1),
                                              oma=c(0,0,1,0),
                                              cex.main=0.8)))


###################################################
### code chunk number 2: Collections.Rnw:84-86
###################################################
data(pHWebs)
pHWebs


###################################################
### code chunk number 3: Collections.Rnw:90-105
###################################################
length(pHWebs)
is.list(pHWebs)
names(pHWebs)

# Access first community in the collection
pHWebs[[1]]

# Access a community by name
pHWebs[['Broadstone']]

# The number of trophic links in Broadstone
NumberOfTrophicLinks(pHWebs[['Broadstone']])

# The number of trophic links in each of the ten webs
sapply(pHWebs, 'NumberOfTrophicLinks')


###################################################
### code chunk number 4: Collections.Rnw:112-114 (eval = FALSE)
###################################################
## length(pHWebs) <- 2 # You can't do this
## pHWebs[1] <- "This will not work"


###################################################
### code chunk number 5: Collections.Rnw:119-120
###################################################
all(FALSE==duplicated(names(pHWebs)))


###################################################
### code chunk number 6: Collections.Rnw:130-138
###################################################
# Returns a new CommunityCollection that contains every other web
pHWebs[seq(1, 10, by=2)]

# Returns a new CommunityCollection with the order reversed
pHWebs[10:1]

# Returns a new CommunityCollection containing only these two webs
pHWebs[c('Old Lodge','Bere Stream')]


###################################################
### code chunk number 7: Collections.Rnw:144-145
###################################################
CollectionCPS(pHWebs)


###################################################
### code chunk number 8: Collections.Rnw:166-168
###################################################
res <- CollectionCPS(pHWebs, properties=c('pH', 'NumberOfNodes'))
res


###################################################
### code chunk number 9: Collections.Rnw:172-174
###################################################
model <- lm(NumberOfNodes ~ pH, data=res)
model


###################################################
### code chunk number 10: Collections.Rnw:177-178
###################################################
summary(model)


###################################################
### code chunk number 11: Collections.Rnw:187-189
###################################################
getOption("SweaveHooks")[["fig"]]()
with(res, plot(pH, NumberOfNodes, pch=19))
abline(model)


###################################################
### code chunk number 12: Collections.Rnw:199-204
###################################################
CollectionCPS(pHWebs, c('pH',
                        'NumberOfNodes',
                        'NumberOfTrophicLinks', 
                        'DirectedConnectance',
                        'NvMSlope'))


###################################################
### code chunk number 13: Collections.Rnw:207-212
###################################################
CollectionCPS(pHWebs, c('pH',
                        S='NumberOfNodes',
                        L='NumberOfTrophicLinks', 
                        C='DirectedConnectance',
                        Slope='NvMSlope'))


###################################################
### code chunk number 14: Collections.Rnw:219-225
###################################################
CollectionCPS(pHWebs, c('pH',
                        S='NumberOfNodes',
                        L='NumberOfTrophicLinks', 
                        C='DirectedConnectance',
                        Slope='NvMSlope',
                        'SumBiomassByClass'))


###################################################
### code chunk number 15: Collections.Rnw:229-235
###################################################
CollectionCPS(pHWebs, c('pH',
                        S='NumberOfNodes',
                        L='NumberOfTrophicLinks', 
                        C='DirectedConnectance',
                        Slope='NvMSlope',
                        B='SumBiomassByClass'))


###################################################
### code chunk number 16: Collections.Rnw:244-250
###################################################
CollectionCPS(pHWebs, list('pH',
                           S='NumberOfNodes',
                           L='NumberOfTrophicLinks', 
                           C='DirectedConnectance',
                           Slope='NvMSlope',
                           B=list('SumBiomassByClass', na.rm=TRUE)))


###################################################
### code chunk number 17: Collections.Rnw:253-257
###################################################
CollectionCPS(pHWebs, c(Basal='FractionBasalNodes', 
                        Intermediate='FractionIntermediateNodes', 
                        TopLevel='FractionTopLevelNodes', 
                        Isolated='FractionIsolatedNodes'))


###################################################
### code chunk number 18: Collections.Rnw:264-272
###################################################
getOption("SweaveHooks")[["fig"]]()
properties <- CollectionCPS(pHWebs, c(S='NumberOfNodes', 
                                      L='NumberOfTrophicLinks', 
                                      'LinkageDensity', 
                                      C='DirectedConnectance'))
par(mfrow=c(1,3))
with(properties, plot(S, L, pch=19))
with(properties, plot(S, LinkageDensity, pch=19))
with(properties, plot(S, C, pch=19))


###################################################
### code chunk number 19: Collections.Rnw:281-282
###################################################
head(CollectionNPS(pHWebs))


###################################################
### code chunk number 20: Collections.Rnw:286-294
###################################################
# A subset of first-class properties
head(CollectionNPS(pHWebs, 'M'))

# Several properties
head(CollectionNPS(pHWebs, c('M','N','Biomass','Degree','IsBasalNode')))

# Named properties
head(CollectionNPS(pHWebs, c('M','N',B='Biomass', 'Degree', Basal='IsBasalNode')))


###################################################
### code chunk number 21: Collections.Rnw:300-301
###################################################
head(CollectionTLPS(pHWebs))


###################################################
### code chunk number 22: Collections.Rnw:304-305
###################################################
head(CollectionTLPS(pHWebs, 'M'))


###################################################
### code chunk number 23: Collections.Rnw:308-309
###################################################
head(CollectionTLPS(pHWebs, c('M','N','Biomass','Degree','IsBasalNode')))


###################################################
### code chunk number 24: Collections.Rnw:312-313
###################################################
head(CollectionTLPS(pHWebs, c('M','N', B='Biomass', D='Degree', Basal='IsBasalNode')))


###################################################
### code chunk number 25: Collections.Rnw:324-325
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pHWebs)


###################################################
### code chunk number 26: Collections.Rnw:331-332
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pHWebs, xlim=c(-14,6), ylim=c(-3,13))


###################################################
### code chunk number 27: Collections.Rnw:340-341
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pHWebs, plot.fn=PlotWebByLevel)


###################################################
### code chunk number 28: Collections.Rnw:346-347
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pHWebs, plot.fn=PlotWebByLevel, ylim=c(1, 4.5))


###################################################
### code chunk number 29: Collections.Rnw:358-360
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pHWebs, plot.fn=PlotNPS, X='Log10M', Y='PreyAveragedTrophicLevel', 
     show.web=FALSE, highlight.nodes=NULL, xlim=c(-14,6), ylim=c(1,4.2))


###################################################
### code chunk number 30: Collections.Rnw:365-367
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pHWebs, plot.fn=PlotTLPS, X='consumer.Log10M', 
     Y='resource.Log10M', xlim=c(-2.5, 5.5), ylim=c(-13.8, 5.5))


###################################################
### code chunk number 31: Collections.Rnw:404-408
###################################################
# Bere Stream has some isolated nodes
CollectionCPS(pHWebs, 'FractionIsolatedNodes')
pHWebs.no.iso <- CollectionApply(pHWebs, RemoveIsolatedNodes)
CollectionCPS(pHWebs.no.iso, 'FractionIsolatedNodes')   # All 0


###################################################
### code chunk number 32: Collections.Rnw:412-416
###################################################
# The number of cannibals in each community
sapply(pHWebs, function(community) length(Cannibals(community)))
pHWebs.no.can <- CollectionApply(pHWebs, RemoveCannibalisticLinks)
sapply(pHWebs.no.can, function(community) length(Cannibals(community)))


###################################################
### code chunk number 33: Collections.Rnw:421-424
###################################################
head(CollectionNPS(pHWebs))
pHWebs.by.M <- CollectionApply(pHWebs, OrderCommunity, 'M')
head(CollectionNPS(pHWebs.by.M))


###################################################
### code chunk number 34: Collections.Rnw:427-429
###################################################
pHWebs.by.M <- CollectionApply(pHWebs, OrderCommunity, 'M', na.last=FALSE)
head(CollectionNPS(pHWebs.by.M))


###################################################
### code chunk number 35: Collections.Rnw:437-439
###################################################
pHWebs.decreasing.pH <- OrderCollection(pHWebs, 'pH', decreasing=TRUE)
CollectionCPS(pHWebs.decreasing.pH)


###################################################
### code chunk number 36: Collections.Rnw:442-444
###################################################
pHWebs.name <- OrderCollection(pHWebs, 'title')
CollectionCPS(pHWebs.name)


###################################################
### code chunk number 37: Collections.Rnw:447-449
###################################################
pHWebs.n.nodes <- OrderCollection(pHWebs, 'NumberOfNodes')
CollectionCPS(pHWebs.n.nodes, c('pH', 'lat', 'NumberOfNodes'))


###################################################
### code chunk number 38: Collections.Rnw:454-456
###################################################
pHWebs.n.nodes.and.lat <- OrderCollection(pHWebs, 'NumberOfNodes', 'lat')
CollectionCPS(pHWebs.n.nodes.and.lat, c('pH', 'lat', 'NumberOfNodes'))


###################################################
### code chunk number 39: Collections.Rnw:468-471
###################################################
data(Millstream)
Millstream
names(Millstream)


###################################################
### code chunk number 40: Collections.Rnw:475-477
###################################################
nps <- CollectionNPS(Millstream)
nps['Synorthocladius sp.'==nps$node,c('community','M','N')]


###################################################
### code chunk number 41: Collections.Rnw:481-485
###################################################
aggregation1 <- AggregateCommunities(Millstream, weight.by='N')

# Satisfy ourselves that each node has been included in the aggregated community
all(sort(unique(nps$node))==sort(NPS(aggregation1)$node))


###################################################
### code chunk number 42: Collections.Rnw:489-490
###################################################
NPS(aggregation1)['Synorthocladius sp.',c('M','N')]


###################################################
### code chunk number 43: Collections.Rnw:493-499
###################################################
# Arithmetic mean of N
mean(nps['Synorthocladius sp.'==nps$node,'N'])

# N-weighted mean of M
weighted.mean(nps['Synorthocladius sp.'==nps$node,'M'], 
              nps['Synorthocladius sp.'==nps$node,'N'])


###################################################
### code chunk number 44: Collections.Rnw:504-513
###################################################
aggregation2 <- AggregateCommunities(Millstream, weight.by=NULL)

NPS(aggregation2)['Synorthocladius sp.',c('M','N')]

# Arithmetic mean of M
mean(nps['Synorthocladius sp.'==nps$node,'M'])

# Arithmetic mean of N
mean(nps['Synorthocladius sp.'==nps$node,'N'])


###################################################
### code chunk number 45: Collections.Rnw:520-523
###################################################
tlps <- CollectionTLPS(Millstream)
tlps['Synorthocladius sp.'==tlps$resource | 
     'Synorthocladius sp.'==tlps$consumer,]


###################################################
### code chunk number 46: Collections.Rnw:526-527
###################################################
TrophicLinksForNodes(aggregation1, 'Synorthocladius sp.')


###################################################
### code chunk number 47: Collections.Rnw:533-536
###################################################
CollectionCPS(Millstream)

data.frame(CPS(aggregation1))


###################################################
### code chunk number 48: Collections.Rnw:545-546
###################################################
CollectionCPS(pHWebs[c('Duddon Pike Beck', 'Mosedal Beck')])


###################################################
### code chunk number 49: Collections.Rnw:550-551
###################################################
CollectionCPS(AggregateCommunitiesBy(pHWebs, 'lat'))


###################################################
### code chunk number 50: Collections.Rnw:566-582
###################################################
data(TL84, TL86)
TL <- CommunityCollection(list(TL84, TL86))
# TL.aggregated is a new Community object containing every species in the TL
all.TL <- AggregateCommunities(TL)    

# Generate a factor of categories
nps <- NPS(all.TL, c('node', 'category'))
categories <- factor(nps$category, levels=c('producer', 'invertebrate', 
                                            'vert.ecto'))

# Order all.TL by categories
all.TL <- OrderCommunity(all.TL, new.order=order(categories, nps$node))

# Create the mapping from node name to ID
map <- 1:NumberOfNodes(all.TL)
names(map) <- unname(NP(all.TL, 'node'))


###################################################
### code chunk number 51: Collections.Rnw:587-591
###################################################
data.frame(ID=1:NumberOfNodes(all.TL), 
           NPS(all.TL, c(Species='node', Category='category', 
                         'M', 'N', TL='PreyAveragedTrophicLevel')), 
           row.names=NULL)


###################################################
### code chunk number 52: Collections.Rnw:603-610
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,2))
for(community in TL)
{
    PlotNvM(community, show.nodes.as='labels', show.web=FALSE, 
            node.labels=map[NP(community, 'node')], xlim=c(-14, 0), 
            ylim=c(-2, 10))
}


