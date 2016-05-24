### R code from vignette source 'kobe.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prelim
###################################################
library(ggplot2)
library(plyr)
library(kobe)
library(reshape)

data(prj1)
data(prj)
data(assmt)

TACs=seq(15000,35000,5000)


###################################################
### code chunk number 2: kobe.Rnw:147-150 (eval = FALSE)
###################################################
## ### Results from ASPIC bootstrapped assessment
## bio   ="http://www.iccat.int/stocka/Models/ASPIC/albs/2011/run2/aspic.bio"
## assmt =kobeAspic(bio)


###################################################
### code chunk number 3: kobe.Rnw:152-153
###################################################
head(assmt)


###################################################
### code chunk number 4: kobe.Rnw:158-161 (eval = FALSE)
###################################################
## ## Results from an ASPIC Projection
## prb  ="http://www.iccat.int/stocka/Models/ASPIC/albs/2011/run2/aspic_15000.prb"
## prj1 =kobeAspic(bio,prb)


###################################################
### code chunk number 5: kobe.Rnw:164-165
###################################################
tail(prj1)


###################################################
### code chunk number 6: kobe.Rnw:179-183 (eval = FALSE)
###################################################
## ## Projections
## TACs=seq(15000,35000,5000)
## prb ="http://www.iccat.int/stocka/Models/ASPIC/albs/2011/run2/aspic_"
## prb =paste(prb,TACs,".prb",sep="")


###################################################
### code chunk number 7: kobe.Rnw:187-189 (eval = FALSE)
###################################################
## ## Results
## prj=kobeAspic(bio,prb,what=c("pts","trks","smry"))


###################################################
### code chunk number 8: kobe.Rnw:192-194
###################################################
class(prj)
names(prj)


###################################################
### code chunk number 9: kobe.Rnw:197-199 (eval = FALSE)
###################################################
## ## add TAC column to data.frame
## prj=llply(prj, transform, TAC=TACs[X1])


###################################################
### code chunk number 10: kobe.Rnw:204-205
###################################################
head(prj$trks)


###################################################
### code chunk number 11: kobe.Rnw:210-211
###################################################
head(prj$smry)


###################################################
### code chunk number 12: kobe.Rnw:216-217
###################################################
data(sims)


###################################################
### code chunk number 13: kobe.Rnw:240-244
###################################################
ggplot(assmt)                                       + 
  geom_hline(aes(yintercept=1),col="red",size=2)    + 
  geom_line( aes(year,stock,group=iter,col=iter))   +
  theme(legend.position="none")


###################################################
### code chunk number 14: kobe.Rnw:258-265
###################################################
### tracks
ggplot(subset(prj$trks,year<=2020))                           +
  geom_line(aes(year,stock,  linetype=Percentile),col="blue") +
  geom_line(aes(year,harvest,linetype=Percentile),col= "red") +
  scale_linetype_manual(values=c(2,1,2))                      +
  coord_cartesian(ylim=c(0,3))                                +
  facet_wrap(~TAC,ncol=2)


###################################################
### code chunk number 15: kobe.Rnw:283-286
###################################################
kp=kobePhase(subset(sims, year==2010 & TAC==15000)) +
         geom_point(aes(stock,harvest,group=Run,col=Run)) 
kp


###################################################
### code chunk number 16: kobe.Rnw:297-299
###################################################
data(sims)
head(sims)


###################################################
### code chunk number 17: kobe.Rnw:311-315
###################################################
dat =subset(sims,year<=2010 & TAC==15000)
trks=ddply(dat,.(Run,year,TAC), function(x) kobeTrks(x$stock,x$harvest,prob=c(0.5)))

head(trks)


###################################################
### code chunk number 18: kobe.Rnw:324-328
###################################################
kp + geom_path( aes(stock,harvest,group=Run,col=Run), data=trks) +
     geom_point(aes(stock,harvest,group=Run), data=subset(trks,year==2010),col="cyan",size=3)+
     facet_wrap(~Run) + 
     theme(legend.position = "none")


###################################################
### code chunk number 19: kobe.Rnw:338-343 (eval = FALSE)
###################################################
## kp2 = kp + geom_path(aes(x,y,group=level),colour="blue",
##                     data=ddply(subset(sims,year==2010 & TAC==15000),.(Run), 
##                                function(pts) kobeProb(pts$stock,pts$harvest,prob=c(0.7,.5,.25)))) +
##                     facet_wrap(~Run) + 
##                     theme(legend.position = "none")


###################################################
### code chunk number 20: kobe.Rnw:347-348 (eval = FALSE)
###################################################
## print(kp2)


###################################################
### code chunk number 21: kobe.Rnw:358-363
###################################################
pts =subset(sims, year==2010 & TAC==15000)

# stock density plot
ggplot(pts) + 
  geom_density(aes(x=stock, y= ..count.., group=Run, fill=Run, alpha=0.4))


###################################################
### code chunk number 22: kobe.Rnw:371-374
###################################################
ggplot(pts) + 
  geom_density(aes(x=stock, y=..count.., group=Run, fill=Run), 
                        fill="grey", col=grey(.9), position = "stack") 


###################################################
### code chunk number 23: kobe.Rnw:381-383
###################################################
### Bespoke Stuff ###
print(kobePhaseMar(transform(pts,group=Run)))          


###################################################
### code chunk number 24: kobe.Rnw:394-408
###################################################
### Pies ###
pie.dat=ddply(subset(sims,year==2010 & TAC==15000),.(Run),kobeSmry,o=T)
pie.dat=ddply(melt(pie.dat,id.vars="Run"),.(Run,variable), 
              function(x) data.frame(value=mean(x$value)))

## pie charts
ggplot(subset(pie.dat,value>0), aes(x =factor(1), y=value, fill = variable)) + 
  geom_bar(width=1,stat="identity") + 
  coord_polar(theta="y") +
  labs(fill='Kobe Quadrant') + xlab('') + ylab('')       +
  scale_fill_manual(values=c("red","green","yellow"))    + 
  facet_wrap(~Run)                                       + 
  scale_x_discrete(breaks=NULL)                          +
  scale_y_continuous(breaks=NULL) 


###################################################
### code chunk number 25: kobe.Rnw:417-457
###################################################
library(akima)
Interp=function(x,levels=seq(0.0,1.0,0.05),
               col   =c(colorRampPalette(c("red4","red"))(12),colorRampPalette(c("yellowgreen","darkgreen"))(8)),
               nIterp=101){

  x=x[!is.na(x[,1]) & !is.na(x[,2]) & !is.na(x[,3]),]
  
  ##### smooth
  t.<-interp(x[,1],x[,2],x[,3],
                    xo=seq(min(x[,1]),   max(x[,1]), length=nIterp),
                    yo=seq(min(x[,2]),   max(x[,2]), length=nIterp),
                    duplicate="mean")
  
  
  res=cbind(expand.grid(x=t.$x,y=t.$y),z=cut(t.$z,levels,include.lowest=T),w=c(t.$z))
  res$col=col[as.numeric(res$z)]
  
  res}

kobe2012=subset(sims,year %in% 2013:2022)
  
pdat=subset(ddply(kobe2012,.(year,TAC),kobeSmry),
            select=c(year,TAC,green,underFished,underFishing))
pdat=melt(pdat,id.vars=c("year","TAC"))
pdat=ddply(pdat, .(variable), function(x) Interp(data.frame(x$year,x$TAC,x$value)))

col.=c(colorRampPalette(c("red4","red"))(12),
       colorRampPalette(c("yellowgreen","darkgreen"))(8))

k2p = ggplot(aes(x=x,y=y,z=w),data=pdat)                 +
           geom_tile(aes(x,y,fill=z))                    +
           scale_fill_manual(values=col.,guide="none")   +
           stat_contour(aes(colour= ..level..),size=1.2,  
                            breaks=c(0.6,0.7,0.8,0.9))   +
           scale_colour_gradient(low="grey", high="black", 
                                 breaks=c(0.6,0.7,0.8,0.9),
                                  labels=c(0.6,0.7,0.8,0.9),limits=c(0.6,1))    +
           facet_wrap(~variable,ncol=1)                       +
           xlab("Year")+ylab("TAC") 
k2p


###################################################
### code chunk number 26: kobe.Rnw:477-483
###################################################
t.=ddply(subset(sims,year %in% 2013:2022),.(year,TAC),  kobeSmry)

k2smTab=list()
k2smTab[[1]]=cast(subset(t., year %in% 2013:2022),TAC~year,value="underFishing")
k2smTab[[2]]=cast(subset(t., year %in% 2013:2022),TAC~year,value="underFished")
k2smTab[[3]]=cast(subset(t., year %in% 2013:2022),TAC~year,value="green")


