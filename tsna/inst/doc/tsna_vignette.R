## ----setup, include=FALSE------------------------------------------------

library(knitr)
knitr::opts_chunk$set(comment='')
knitr::opts_chunk$set(cache=FALSE)

## ------------------------------------------------------------------------
library(tsna)

## ------------------------------------------------------------------------
?tsna
?tPath

## ------------------------------------------------------------------------
library(networkDynamicData)
library(sna)

## ------------------------------------------------------------------------
data(moodyContactSim)
plot(moodyContactSim,displaylabels = TRUE,main='aggregate network')

## ------------------------------------------------------------------------
as.data.frame(moodyContactSim)

## ------------------------------------------------------------------------
coords<-plot(moodyContactSim,
     displaylabels=TRUE,
     label.cex=0.8,
     label.pos=5,
     vertex.col='white',
     vertex.cex=3,
     edge.label=sapply(get.edge.activity(moodyContactSim),function(e){
       paste('(',e[,1],'-',e[,2],')',sep='')
     }),
     edge.label.col='blue',
     edge.label.cex=0.7
   )

## ------------------------------------------------------------------------
par(mfcol=c(1,2))
plot(network.extract(moodyContactSim,at=215),
     main='network at time 215',
     displaylabels=TRUE,
     label.cex=0.6,
     label.pos=5,
     vertex.col='white',
     vertex.cex=3,
     coord=coords)
plot(network.extract(moodyContactSim,at=750),
     main='network at time 750',
     displaylabels=TRUE,
     label.cex=0.6,
     label.pos=5,
     vertex.col='white',
     vertex.cex=3,
     coord=coords)
par(mfcol=c(1,1))


## ------------------------------------------------------------------------
v10path<-tPath(moodyContactSim,v=10,start=0)
print(v10path)

## ------------------------------------------------------------------------
plot(v10path,coord=coords, displaylabels=TRUE)

## ------------------------------------------------------------------------
v1path<-tPath(moodyContactSim,v=1,start=0)

## ------------------------------------------------------------------------
par(mfcol=c(1,2)) # set up side-by-side plot

plotPaths(moodyContactSim,v10path,coord=coords, main='fwd path from v10')
plotPaths(moodyContactSim,v1path,coord=coords, main = 'fwd path from v1')

par(mfcol=c(1,1)) # turn off side-by-side plots

## ------------------------------------------------------------------------
# or draw both plots on the  them both on the same network
plotPaths(moodyContactSim,coord=coords,list(v10path,v1path))

## ------------------------------------------------------------------------
plotPaths(moodyContactSim, list(
          tPath(moodyContactSim,v=10,direction='bkwd',type='latest.depart'),
          tPath(moodyContactSim,v=10)))


## ------------------------------------------------------------------------
par(mfcol=1:2)
plot(tPath(moodyContactSim,v=1,start=0),coord=coords,
     main='tPath from v1 @ t=0')
plot(tPath(moodyContactSim,v=1,start=500),coord=coords,
     main='tPath from v1 @ t=500')
par(mfcol=c(1,1))

## ------------------------------------------------------------------------
pathCompare<-network.initialize(7)
network.vertex.names(pathCompare)<-LETTERS[1:7]
add.edges.active(pathCompare,tail=c(1,2),head=c(2,7),onset=c(1,4),terminus=c(2,5))
add.edges.active(pathCompare,tail=c(1,3),head=c(3,7),onset=c(0,6),terminus=c(2,7))
add.edges.active(pathCompare,tail=c(1,4),head=c(4,7),onset=c(4,5),terminus=c(5,6))
add.edges.active(pathCompare,tail=c(1,5),head=c(5,7),onset=c(6,9),terminus=c(7,10))
add.edges.active(pathCompare,tail=c(1,6),head=c(6,7),onset=c(4,10),terminus=c(5,11))
as.data.frame(pathCompare)

## ----fig.width=8,fig.height=8--------------------------------------------
# pre-define some coords for arbitrary positioning
coords<-cbind(c(0,0.5,0.5,0.5,0.5,0.5,1),c(0.3,0.15,0.3,0.45,0.65,0.8,0.7))
# do the plot
plot(pathCompare,
     coord=coords,jitter=FALSE,
     #mode='circle',
     displaylabels=TRUE, vertex.col='white',
     edge.label=get.edge.activity(pathCompare),edge.label.cex=0.8,
     edge.lwd=4,
     edge.col=c('blue','blue','red','red','green','green','orange','orange','pink','pink'),
    main='Comparison of fwd temporal path types from A to G')
# plot a legend
legend(-0.3,1,legend = c('earliest leaving (ACG @ 6)',
                         'earliest arriving (ABG @ 4)',
                         'latest leaving (AEG @ 10)',
                         'quickest (ADG @ 5)',
                         'latest arriving (AFG @ 11)'),
               fill=c('red','blue','orange','green','pink'),
       cex=0.8)

## ------------------------------------------------------------------------
library(networkDynamicData)
data(concurrencyComparisonNets)

## ------------------------------------------------------------------------
baseTrees<-tReach(base,sample=25)
monogTrees<-tReach(monog,sample=25)
middleTrees<-tReach(middle,sample=25)

## ------------------------------------------------------------------------
baseTrees

## ------------------------------------------------------------------------
monogTrees

## ------------------------------------------------------------------------
boxplot(cbind(baseTrees,monogTrees,middleTrees),
        main='fwd-reachable set size distributions for nets of varying concurrency')

## ------------------------------------------------------------------------
hist(baseTrees, main='fwd-reach size distributions',
     ylim=c(0,50),xlim=c(0,1000),
     breaks=seq(from=0,to=1000,by=50),
     col='#55000033',xlab='reachable set size')
hist(monogTrees,ylim=c(0,50),xlim=c(0,1000),
     breaks=seq(from=0,to=1000,by=50),
     col='#00550033',add=TRUE)
hist(middleTrees,ylim=c(0,50),xlim=c(0,1000),
     breaks=seq(from=0,to=1000,by=50),
     col='#00005533',add=TRUE)
legend(800,50,legend=c('base','monog','middle'),
       fill=c('#55000033','#00550033','#00005533'))

## ------------------------------------------------------------------------
mean(degree(as.network(base)))
mean(degree(as.network(monog)))
mean(degree(as.network(middle)))

## ----fig.width=8,fig.height=8--------------------------------------------
nets4<-replicate(4,list(network(matrix(rbinom(16,5,0.1),ncol=4,nrow=4))))
par(mfcol=c(2,2))
plot(nets4[[1]],displaylabels=TRUE,main='t0')
plot(nets4[[2]],displaylabels=TRUE,main='t1')
plot(nets4[[3]],displaylabels=TRUE,main='t2')
plot(nets4[[4]],displaylabels=TRUE,main='t3')
par(mfcol=c(1,1))

## ------------------------------------------------------------------------
nets4Dyn<-networkDynamic(network.list=nets4)
nets4Projected<-timeProjectedNetwork(nets4Dyn)

## ------------------------------------------------------------------------
network.size(nets4Projected)
network.vertex.names(nets4Projected)

## ----fig.width=8,fig.height=8--------------------------------------------
plot(nets4Projected,
     displaylabels=TRUE,
     mode='kamadakawai',
     edge.col=ifelse(nets4Projected%e%'edge.type'=='identity_arc','gray','black'))

## ----fig.width=8,fig.height=8--------------------------------------------
changes<-get.change.times(moodyContactSim)
moodyProj<-timeProjectedNetwork(moodyContactSim,onsets=changes,termini=changes)
plot(moodyProj,
     mode='kamadakawai',
     vertex.cex=0.3,
     arrowhead.cex=0.1,
     edge.col=ifelse(moodyProj%e%'edge.type'=='identity_arc','gray','black'))

## ------------------------------------------------------------------------
plot(tEdgeDissolution(base),main="Edge dissolution counts for network 'base'")
plot(tEdgeFormation(base), main="Edge formation counts for network 'base'")

## ------------------------------------------------------------------------
tEdgeDissolution(base,result.type = 'fraction',time.interval = 10)

## ------------------------------------------------------------------------
tEdgeFormation(base,result.type = 'fraction',time.interval = 10)

## ------------------------------------------------------------------------
data(harry_potter_support)

## ------------------------------------------------------------------------
tSnaStats(harry_potter_support,snafun='gtrans')

## ------------------------------------------------------------------------
# compute triad census scores for each time point
tSnaStats(harry_potter_support,snafun='triad.census')

## ------------------------------------------------------------------------
# compute degrees
bet<-tSnaStats(harry_potter_support,snafun='betweenness')
nrow(bet)
ncol(bet)
bet[,25,drop=FALSE]
class(bet)

## ------------------------------------------------------------------------
colMeans(bet,na.rm = TRUE)

## ------------------------------------------------------------------------
rowMeans(bet)

## ------------------------------------------------------------------------
prestScores<-tSnaStats(base,'prestige',time.interval=25,rescale=TRUE)

## ------------------------------------------------------------------------
 tErgmStats(base,'~edges+concurrent',
               start=0,end=100,time.interval = 10)

## ------------------------------------------------------------------------
 plot(
   tErgmStats(base,'~edges+concurrent',
                start=0,end=100,time.interval = 10)
    )

## ------------------------------------------------------------------------
edgeDuration(moodyContactSim)

## ------------------------------------------------------------------------
summary(edgeDuration(moodyContactSim))
hist(edgeDuration(moodyContactSim))

## ------------------------------------------------------------------------
data(concurrencyComparisonNets)
hist(edgeDuration(base),ylim=c(0,800))
hist(edgeDuration(middle),ylim=c(0,800))
hist(edgeDuration(monog),ylim=c(0,800))

## ------------------------------------------------------------------------
which(edgeDuration(monog,mode='counts')>1)
which(edgeDuration(moodyContactSim,mode='counts')>1)

## ------------------------------------------------------------------------
edgeDuration(monog,e=valid.eids(monog)[105],subject='edges')
edgeDuration(monog,e=valid.eids(monog)[105],subject='spells')

## ------------------------------------------------------------------------
mean(edgeDuration(base,subject = 'edges'))
mean(edgeDuration(base,subject = 'spells'))

## ------------------------------------------------------------------------
data(windsurfers)
vertexDuration(windsurfers)
table(vertexDuration(windsurfers))
hist(vertexDuration(windsurfers))
hist(vertexDuration(windsurfers,subject='spells'))

## ------------------------------------------------------------------------
data(McFarland_cls33_10_16_96)
tiedDuration(cls33_10_16_96, mode='counts')

## ------------------------------------------------------------------------
cls33_10_16_96%v%'type'

## ------------------------------------------------------------------------
tiedDuration(cls33_10_16_96, mode='counts',neighborhood = 'in')

## ------------------------------------------------------------------------
plot(tiedDuration(cls33_10_16_96, mode='counts',neighborhood = 'out'),
     tiedDuration(cls33_10_16_96, mode='counts',neighborhood = 'in'),
     xlab='# speaking events',ylab='# spoken to events',main='McFarland classroom network, speaking vs. spoken to' )
text(tiedDuration(cls33_10_16_96, mode='counts',neighborhood = 'out'),
     tiedDuration(cls33_10_16_96, mode='counts',neighborhood = 'in'),
     label=cls33_10_16_96%v%'type',cex=0.8,pos=4)

## ------------------------------------------------------------------------
plot(sort(tiedDuration(base)),type='l',ylim=c(0,400),
     main='sorted tiedDuration for concurrency scenearios',
     xlab='sorted vertices',ylab='duration that each vertex is connected', col ='#55000033',lwd=4)
points(sort(tiedDuration(monog)),type='l',col='#00550033',lwd=4)
points(sort(tiedDuration(middle)),type='l',col='#00005533',lwd=4)
legend(200,300,legend=c('base','monog','middle'),
       fill=c('#55000033','#00550033','#00005533'))

mean(tiedDuration(base))
mean(tiedDuration(monog))
mean(tiedDuration(middle))

## ----fig.width=8,fig.height=4--------------------------------------------
par(mfcol=c(1,3))
plot(degree(as.network(base)),tiedDuration(base),xlim=c(0,25),ylim=c(0,350),main='base')
plot(degree(as.network(middle)),tiedDuration(middle),xlim=c(0,25), ylim=c(0,350),main='middle')
plot(degree(as.network(monog)),tiedDuration(monog),xlim=c(0,25),ylim=c(0,350),main='monog')
par(mfcol=c(1,1))

## ------------------------------------------------------------------------
data(hospital_contact)
plot(degree(as.network(hospital),gmode = 'graph'),tiedDuration(hospital),
     xlab='aggregate degree (total number of unique contacts)',
     ylab='total contact duration (seconds)',
     main='Vertices in hospital RFID proximity contact network')

## ------------------------------------------------------------------------
data(moodyContactSim)
data(harry_potter_support)
data(onlineNetwork)
data(vanDeBunt_students)
data(McFarland_cls33_10_16_96)
data(windsurfers)
data(hospital_contact)
data(concurrencyComparisonNets)
nets<-list(
  moodyContactSim=moodyContactSim,
  hospital=hospital,
  base=base,
  monog=monog,
  harry_potter=harry_potter_support,
  onlineNet=onlineNet,
  vanDeBunt=vanDeBunt_students,
  McFarland=cls33_10_16_96,
  windsurfers=windsurfers)

## ------------------------------------------------------------------------
par(mfcol=c(3,3))
for (n in seq_along(nets)){
  hist(edgeDuration(nets[[n]]),main=names(nets)[n],xlab='duration')
}
par(mfcol=c(1,1))

## ------------------------------------------------------------------------
par(mfcol=c(3,3))
for (n in seq_along(nets)){
  hist(edgeDuration(nets[[n]],mode = 'counts'),main=names(nets)[n],xlab='duration')
}
par(mfcol=c(1,1))

## ------------------------------------------------------------------------
tEdgeDensity(base)

## ------------------------------------------------------------------------
tEdgeDensity(base,agg.unit = 'dyad')

## ------------------------------------------------------------------------
edd<-sapply(nets,tEdgeDensity)
plot(edd,main='edge duration density',xaxt='n',xlab='networks')
text(edd,label=names(edd),pos=4)

## ------------------------------------------------------------------------
eed<-sapply(nets,tEdgeDensity,mode='event')
plot(eed,main='edge event density',xaxt='n',xlab='networks')
text(eed,label=names(eed),pos=4)

## ------------------------------------------------------------------------
ddd<-sapply(nets,tEdgeDensity,agg.unit='dyad')
plot(ddd,main='dyad duration density',xaxt='n',xlab='networks')
text(ddd,label=names(ddd),pos=4)

## ------------------------------------------------------------------------
data(McFarland_cls33_10_16_96)
pShiftCount(cls33_10_16_96)

## ------------------------------------------------------------------------
sliceCounts<- lapply(seq(from = 0,to=45,by = 5),function(onset){
  pShiftCount(network.extract(cls33_10_16_96,onset,length = 5))
})
# convert to a matrix
sliceCounts<-do.call(rbind,sliceCounts)
sliceCounts
# plot
plot(sliceCounts[,'AB-BA'],type='b',col='blue',ylim=c(0,170),
     main='pShift counts for 5-min intervals of cls33',
     ylab='count of selected pShift',xlab='slice index')
points(sliceCounts[,'AB-AY'],type='b',col='red')
points(sliceCounts[,'AB-XA'],type='b',col='green')

## ------------------------------------------------------------------------
citation('tsna')

