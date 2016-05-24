### R code from vignette source 'usageExamples.Rnw'

###################################################
### code chunk number 1: package
###################################################
options(keep.source = TRUE, width = 60)
packageInfo <- packageDescription("vipor")
library(vipor)
packageKeywords<-"visualization, display, one dimensional, grouped, groups, violin, scatter, points, quasirandom, beeswarm, van der Corput"


###################################################
### code chunk number 2: vpPlot (eval = FALSE)
###################################################
##   library(vipor)
##   set.seed(12345)
##   n<-100
##   dat<-rnorm(n*2)
##   labs<-rep(c('a','b'),n)
##   vpPlot(labs,dat)


###################################################
### code chunk number 3: showVpPlot
###################################################
  library(vipor)
  set.seed(12345)
  n<-100
  dat<-rnorm(n*2)
  labs<-rep(c('a','b'),n)
  vpPlot(labs,dat)


###################################################
### code chunk number 4: vpOpts (eval = FALSE)
###################################################
##   vpPlot(labs,dat,las=1,ylab='Data',col=rep(1:2,n))
##   abline(h=0,lty=2)


###################################################
### code chunk number 5: showVpOpts
###################################################
  vpPlot(labs,dat,las=1,ylab='Data',col=rep(1:2,n))
  abline(h=0,lty=2)


###################################################
### code chunk number 6: vpFactors (eval = FALSE)
###################################################
##   labs2<-factor(labs,levels=c('b','a'))
##   vpPlot(labs2,dat,las=1,ylab='Data',col=rep(1:2,n))
##   abline(h=0,lty=2)


###################################################
### code chunk number 7: showVpFactors
###################################################
  labs2<-factor(labs,levels=c('b','a'))
  vpPlot(labs2,dat,las=1,ylab='Data',col=rep(1:2,n))
  abline(h=0,lty=2)


###################################################
### code chunk number 8: offsetX
###################################################
  offsets<-offsetX(dat,labs)
  head(offsets,4)
  xPos<-vpPlot(labs,dat)
  head(xPos,4)
  xPos2<-rep(1:2,n)+offsets
  head(xPos2,4)
  all(xPos==xPos2)


###################################################
### code chunk number 9: distAdjust (eval = FALSE)
###################################################
##   dat <- list(
##     'Normal'=rnorm(50),
##     'Dense normal'= rnorm(500),
##     'Bimodal'=c(rnorm(100), rnorm(100,5)),
##     'Extremes'=rcauchy(100)
##   )
##   par(mfrow=c(4,1), mar=c(2.5,3.1, 1.2, 0.5),mgp=c(2.1,.75,0),
##   cex.axis=1.2,cex.lab=1.2,cex.main=1.2)
##   dummy<-sapply(names(dat),function(label) {
##     y<-dat[[label]]
##     offsets <- list(
##       'defaults'=offsetX(y),  # Default
##       'adjust=2'=offsetX(y, adjust=2),    # More smoothing
##       'adjust=.1'=offsetX(y, adjust=0.1),  # Tighter fit
##       'width=.1'=offsetX(y, width=0.1),    # Less wide
##       'nbins=100'=offsetX(y, nbins=100)    # Less bins
##     )  
##     ids <- rep(1:length(offsets), each=length(y))
##     plot(unlist(offsets) + ids, rep(y, length(offsets)), ylab='y value',
##       xlab='', xaxt='n', pch=21,
##       col='#00000099',bg='#00000033',las=1,main=label)
##     axis(1, 1:length(offsets), names(offsets))
##   })


###################################################
### code chunk number 10: showDistAdjust
###################################################
  dat <- list(
    'Normal'=rnorm(50),
    'Dense normal'= rnorm(500),
    'Bimodal'=c(rnorm(100), rnorm(100,5)),
    'Extremes'=rcauchy(100)
  )
  par(mfrow=c(4,1), mar=c(2.5,3.1, 1.2, 0.5),mgp=c(2.1,.75,0),
  cex.axis=1.2,cex.lab=1.2,cex.main=1.2)
  dummy<-sapply(names(dat),function(label) {
    y<-dat[[label]]
    offsets <- list(
      'defaults'=offsetX(y),  # Default
      'adjust=2'=offsetX(y, adjust=2),    # More smoothing
      'adjust=.1'=offsetX(y, adjust=0.1),  # Tighter fit
      'width=.1'=offsetX(y, width=0.1),    # Less wide
      'nbins=100'=offsetX(y, nbins=100)    # Less bins
    )  
    ids <- rep(1:length(offsets), each=length(y))
    plot(unlist(offsets) + ids, rep(y, length(offsets)), ylab='y value',
      xlab='', xaxt='n', pch=21,
      col='#00000099',bg='#00000033',las=1,main=label)
    axis(1, 1:length(offsets), names(offsets))
  })


###################################################
### code chunk number 11: varwidth (eval = FALSE)
###################################################
##   dat <- list(
##     '10 points'=rnorm(10),
##     '50 points'=rnorm(50,2),
##     '200 points'=c(rnorm(400), rnorm(100,5)),
##     '5000 points'= rnorm(5000,1)
##   )
##   labs<-rep(names(dat),sapply(dat,length))
##   labs<-factor(labs,levels=unique(labs))
##   vpPlot( labs,unlist(dat),offsetXArgs=list(varwidth=TRUE),
##     las=1,ylab='Value',col='#00000066',bg='#00000022',pch=21)


###################################################
### code chunk number 12: showVarwidth
###################################################
  dat <- list(
    '10 points'=rnorm(10),
    '50 points'=rnorm(50,2),
    '200 points'=c(rnorm(400), rnorm(100,5)),
    '5000 points'= rnorm(5000,1)
  )
  labs<-rep(names(dat),sapply(dat,length))
  labs<-factor(labs,levels=unique(labs))
  vpPlot( labs,unlist(dat),offsetXArgs=list(varwidth=TRUE),
    las=1,ylab='Value',col='#00000066',bg='#00000022',pch=21)


###################################################
### code chunk number 13: vpBeaver (eval = FALSE)
###################################################
##   y<-c(beaver1$temp,beaver2$temp)
##   x<-rep(
##     c('Beaver 1','Beaver 2'),
##     c(nrow(beaver1),nrow(beaver2))
##   )
##   vpPlot(x,y,las=1, ylab='Body temperature',
##     pch=21, col='#00000099',bg='#00000033')


###################################################
### code chunk number 14: showBeaver
###################################################
  y<-c(beaver1$temp,beaver2$temp)
  x<-rep(
    c('Beaver 1','Beaver 2'),
    c(nrow(beaver1),nrow(beaver2))
  )
  vpPlot(x,y,las=1, ylab='Body temperature',
    pch=21, col='#00000099',bg='#00000033')


###################################################
### code chunk number 15: vpGene (eval = FALSE)
###################################################
##   ints<-integrations[integrations$nearestGene>0,]
##   y<-log(ints$nearestGene)
##   x<-as.factor(paste(ints$study,ints$latent))
##   activeCols<-c('Expressed'='#FF000033','Unexpressed'='#0000FF33')
##   cols<-activeCols[ints$latent]
##   par(mar=c(4,7,.1,.1))
##   vpPlot(x,y,las=2, ylab='Log distance to gene',xaxt='n',
##     pch=21, col=cols,bg=cols,cex=.7)
##   uniqX<-levels(x)
##   prettyX<-tapply(1:length(uniqX),sub('(Une|E)xpressed$','',uniqX),mean)
##   axis(1,prettyX,names(prettyX),las=2)
##   legend(grconvertX(0.01,from='ndc'),grconvertY(0.15,from='ndc'),
##     names(activeCols),pch=21,col=cols,pt.bg=activeCols,xpd=NA)


###################################################
### code chunk number 16: showGene
###################################################
  ints<-integrations[integrations$nearestGene>0,]
  y<-log(ints$nearestGene)
  x<-as.factor(paste(ints$study,ints$latent))
  activeCols<-c('Expressed'='#FF000033','Unexpressed'='#0000FF33')
  cols<-activeCols[ints$latent]
  par(mar=c(4,7,.1,.1))
  vpPlot(x,y,las=2, ylab='Log distance to gene',xaxt='n',
    pch=21, col=cols,bg=cols,cex=.7)
  uniqX<-levels(x)
  prettyX<-tapply(1:length(uniqX),sub('(Une|E)xpressed$','',uniqX),mean)
  axis(1,prettyX,names(prettyX),las=2)
  legend(grconvertX(0.01,from='ndc'),grconvertY(0.15,from='ndc'),
    names(activeCols),pch=21,col=cols,pt.bg=activeCols,xpd=NA)


