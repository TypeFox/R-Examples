# The following code should generate all of the figures
# (with the exception of figure 3) that are included 
# in Handcock & Aldrich (2002).  The full report can be 
# obtained at http://www.csss.washington.edu/Papers

library(reldist)
data(nls, package="reldist")

postscript("Fig1.ps", width=4.0,height=7.5,horizontal=FALSE)
par(mfrow=c(2,1))
par(err=-1)
#par(mar=c(2.5,4.1,2.1,2.1), mgp=c(1.5,1,0))
#
# Set the random seed
#
# First original from Splus
#
#.Random.seed <- c(57, 8, 4, 7, 43, 3, 0, 9, 24, 13, 53, 3)
#
set.seed(sum(c(57, 8, 4, 7, 43, 3, 0, 9, 24, 13, 53, 3)))
#
schpermwage1<-sample(original$chpermwage,size=100000,
                     prob=original$wgt/sum(original$wgt),replace=TRUE)
schpermwage2<-sample(recent$chpermwage,size=100000,
                     prob=recent$wgt/sum(recent$wgt),replace=TRUE)
#
# wage change (34 wage - 16 wage)
#
# Try to overlapping densities
#
nbar <- log(length(schpermwage1), base = 2) + 1
kwidth <- diff(range(schpermwage1))/nbar * 0.5
swidth <- min(sqrt(var(schpermwage1)),diff(quantile(schpermwage1,c(0.25,0.75)))/1.34)
swidth <- 1.06*swidth*length(schpermwage1)^(-1/5)
#
# Following is the width I used
# kwidth <- 1.3*kwidth
#
# Use Silverman's width
#
kwidth <- swidth
kwidth <- 0.2597
dens1 <- density(schpermwage1, n = 500, width=kwidth)
plot(x = (dens1$x), y = dens1$y, type = "l", 
  xlab = "change in log permanent wage", ylab = "density", 
  axes = FALSE,
  xlim = c(-1, 3),
  ylim=c(0,1.2))
title(main="(a)",cex=0.6)
axis(side = 1)
axis(side = 2)
fig1legend <- list(x=c(1.2,1.2),y=c(1.2,1.2))                  
legend(fig1legend,lty=1:2,cex=0.5, bty="n",
legend=c("original cohort","recent cohort"))       
#
swidth <- min(sqrt(var(schpermwage2)),diff(quantile(schpermwage2,c(0.25,0.75)))/1.34)
swidth <- 1.06*swidth*length(schpermwage2)^(-1/5)
#
# Following is the width I used
# kwidth <- 1.3*kwidth
#
# Use Silverman's width
#
kwidth <- swidth
kwidth <- 0.455
#
dens2 <- density(schpermwage2, n = 500, width=kwidth)
lines(x = (dens2$x), y = dens2$y, type = "l",lty=2) 
#
# Now do the Lorenz
#
swage1 <- sort(recent$chpermwage)
swage2 <- sort(original$chpermwage)
xout <- (0:1000)/1000
#
alpha <- seq(along=swage2)/length(swage2)
galpha <- cumsum(swage2)/sum(swage2)
fn1 <- approx(x=alpha,y=galpha,xout=xout)
plot(x = alpha, y = galpha, type = "l", 
  xlab = "proportion of population", 
  ylab = "proportion of wages", 
  ylim=c(0,1.0))
legend(x=c(0,0),y=c(1.03,1.03),lty=1:2,cex=0.5, bty="n",
 legend=c("original cohort","recent cohort"))
#
alpha <- seq(along=swage1)/length(swage1)
galpha <- cumsum(swage1)/sum(swage1)
lines(x = alpha, y = galpha, lty = 2)
lines(alpha,alpha,lty=3)
title(main="(b)",cex=0.6)
#
# Now do the RD with smooth=0.4
#
postscript("Fig2.ps", width=4.0,height=7.5,horizontal=FALSE)
par(err=-1)
par(mfrow=c(2,1))
#par(mar=c(2.5,4.1,2.1,2.1), mgp=c(1.5,1,0))
#
fig2a <- reldist(y=recent$chpermwage,yo=original$chpermwage,ci=FALSE,
        smooth=0.4,
	yowgt=original$wgt,ywgt=recent$wgt,
	cdfplot=TRUE,
	yolabs=seq(-1,3,by=0.5),
	ylabs=seq(-1,3,by=0.5),
	cex=0.8,
        ylab="proportion of the recent cohort",
        xlab="proportion of the original cohort")
title(main="(a)",cex=0.6)
#
fig2b <- reldist(y=recent$chpermwage,yo=original$chpermwage,ci=FALSE,
        smooth=0.4,
	yowgt=original$wgt,ywgt=recent$wgt,bar=TRUE,
	yolabs=seq(-1,3,by=0.5),
	ylim=c(0,2.5),cex=0.8,
        ylab="Relative Density",
        xlab="Proportion of the Original Cohort")
title(main="(b)",cex=0.6)
#
# Now do the RD with smooth=1.2
#
postscript("Fig3.ps", width=4.0,height=7.5,horizontal=FALSE)
par(err=-1)
par(mfrow=c(2,1))
#
fig3a <- reldist(y=recent$chpermwage,yo=original$chpermwage,ci=FALSE,
        smooth=1.2, method="loglik",
	yowgt=original$wgt,ywgt=recent$wgt,
	cdfplot=TRUE,
	yolabs=seq(-1,3,by=0.5),
	ylabs=seq(-1,3,by=0.5),
	cex=0.8,
        ylab="proportion of the recent cohort",
        xlab="proportion of the original cohort")
title(main="(a)",cex=0.6)
#
fig3b <- reldist(y=recent$chpermwage,yo=original$chpermwage,ci=FALSE,
        smooth=1.2,
	yowgt=original$wgt,ywgt=recent$wgt,bar=TRUE,
	yolabs=seq(-1,3,by=0.5),
	ylim=c(0,2.5),cex=0.8,
        ylab="Relative Density",
        xlab="Proportion of the Original Cohort")
title(main="(b)",cex=0.6)
#
postscript(file = "Fig4.ps",width=6.5,height=3.0,horizontal=TRUE)
par(err=-1)
par(mfrow=c(1,3))
#
# Overall r.d. of Y to Yo
# This is a measure of the total differences between Y and Yo.
#
# original$chpermwage is the 
# Predicted change in real log wage from age 16 to 34 for Cohort 1 

g10 <- reldist(y=recent$chpermwage, yo=original$chpermwage, smooth=0.4, ci=FALSE,
	ywgt=recent$wgt, yowgt=original$wgt,
        method="gam",
	yolabs=seq(-1,3,by=0.5),
	ylim=c(0.5,3.0),
	bar=TRUE, quiet=FALSE,
	xlab="proportion of the original cohort")
title(main=paste("(a) entropy = ",format(g10$entropy,digits=3)),cex=0.6)
abline(h=1,lty=2)
#
# calculate the r.d. of Y HAD the 
# only difference from Yo been and ADDITIVE location shift.
# This is a measure of the differences between Y and Yo
# DUE to location differences.
#
# additive location effect
#
# RD of Y to YA
#
g1A <- reldist(y=recent$chpermwage, yo=original$chpermwage,
	ywgt=recent$wgt, yowgt=original$wgt,
	show="effect",
	bar=TRUE, quiet=FALSE,
	ylim=c(0.5,3.0), ylab="",
	smooth=0.4, ci=FALSE,
	yolabs=seq(-1,3,by=0.5),
              xlab="proportion of the original cohort")
title(main=paste("(b) entropy = ",format(entropy(g1A,g10),digits=3)),cex=0.6)
abline(h=1,lty=2)
#
# calculate the r.d. of Y with the same ADDITIVE location as Yo
# to Yo.
# This is a measure of the differences between Y and Yo
# NOT DUE to the ADDITIVE difference in location.
#
# additive shape effect
#
# RD of YA to Yo
#
gA0 <- reldist(y=recent$chpermwage, yo=original$chpermwage, smooth=0.4, ci=FALSE,
	ywgt=recent$wgt, yowgt=original$wgt,
	show="residual",
	bar=TRUE, quiet=FALSE,
	ylim=c(0.5,3.0), ylab="",
	yolabs=seq(-1,3,by=0.5),
              xlab="proportion of the original cohort")
title(main=paste("(c) entropy = ",format(gA0$entropy,digits=3)),cex=0.6)
abline(h=1,lty=2)
#
# print out polarizations
#
format(rpy(y=recent$chpermwage,yo=original$chpermwage,
       ywgt=recent$wgt,yowgt=original$wgt,pvalue=TRUE),
       digits=3)
format(rpluy(y=recent$chpermwage,yo=original$chpermwage,
       ywgt=recent$wgt,yowgt=original$wgt,pvalue=TRUE),
       digits=3)
format(rpluy(y=recent$chpermwage,yo=original$chpermwage,
       ywgt=recent$wgt,yowgt=original$wgt,upper=TRUE,pvalue=TRUE),
       digits=3)
#
# endeduc is the final education 
#
e1 <- original$endeduc
e1[e1 < 8] <- 8
e1[e1 > 18] <- 18
e2 <- recent$endeduc
e2[e2 < 8] <- 8
e2[e2 > 18] <- 18
postscript("Fig5.ps", width=4.5,height=4.5,horizontal=FALSE)
par(err=-1)
g10 <- rddist(y=e2, yo=e1, pool=1, ci=FALSE, quiet=FALSE,
        ywgt=recent$wgt,yowgt=original$wgt,
              yolabs=sort(unique(e1)),
              ylab="relative density",
              xlab="proportion of the original cohort")
title(sub=paste("entropy = ",format(entropy(g10),digits=3)))
abline(h=1,lty=2)
#
postscript(file = "Fig6.ps",width=6.5,height=3.0,horizontal=TRUE)
par(err=-1)
par(mfrow=c(1,3))
#
# Overall r.d. of Y to Yo
# This is a measure of the total differences between Y and Yo.
#
# original$chpermwage is the 
# Predicted change in real log wage from age 16 to 34 for Cohort 1 
#
# generate synthetic Recent data in the same age proportions as
# the Original starting ages
#
# Now the Recent compositionally adjusted for attrition
# and age relative to the Original
#
i3x <- sample(seq(along=original$chpermwage), 
        size = 10*length(original$chpermwage), 
        prob=rdsamp(e2,e1,recent$wgt,original$wgt),
		 replace = TRUE)
schpermwage1 <- original$chpermwage[i3x]
wschpermwage1 <- original$wgt[i3x]
#
g10 <- reldist(y=recent$chpermwage, yo=original$chpermwage, smooth=0.4, ci=FALSE,
	ywgt=recent$wgt, yowgt=original$wgt,
	yolabs=seq(-1,3,by=0.5),
	ylim=c(0.5,3.0),
	bar=TRUE, quiet=FALSE,
	xlab="proportion of the original cohort")
title(main=paste("(a) entropy = ",format(g10$entropy,digits=3)),cex=0.6)
abline(h=1,lty=2)
#
# calculate the r.d. of Y HAD the 
# only difference from Yo been and ADDITIVE location shift.
# This is a measure of the differences between Y and Yo
# DUE to location differences.
#
# additive location effect
#
# RD of Y to YA
#
#par(xaxs="d", yaxs="d")
g1A <- reldist(y=schpermwage1, yo=original$chpermwage,
	yowgt=original$wgt, ywgt=wschpermwage1,
	bar=TRUE, quiet=FALSE,
	ylim=c(0.5,3.0), ylab="",
	smooth=0.4, ci=FALSE,
	yolabs=seq(-1,3,by=0.5),
              xlab="proportion of the original cohort")
title(main=paste("(b) entropy = ",format(entropy(g1A,g10),digits=3)),cex=0.6)
abline(h=1,lty=2)
#
# calculate the r.d. of Y with the same ADDITIVE location as Yo
# to Yo.
# This is a measure of the differences between Y and Yo
# NOT DUE to the ADDITIVE difference in location.
#
# additive shape effect
#
# RD of YA to Yo
#
gA0 <- reldist(y=recent$chpermwage, yo=schpermwage1, smooth=0.4, ci=FALSE,
	ywgt=recent$wgt, yowgt=wschpermwage1,
	bar=TRUE, quiet=FALSE,
	ylim=c(0.5,3.0), ylab="",
	yolabs=seq(-1,3,by=0.5),
              xlab="proportion of the original cohort")
title(main=paste("(c) entropy = ",format(gA0$entropy,digits=3)),cex=0.6)
abline(h=1,lty=2)
#
# print out polarizations
#
format(rpy(y=recent$chpermwage,yo=original$chpermwage,
           ywgt=recent$wgt,yowgt=original$wgt,pvalue=TRUE),
       digits=3)
format(rpluy(y=schpermwage1,yo=original$chpermwage,
             ywgt=wschpermwage1,yowgt=original$wgt,
             pvalue=TRUE),
       digits=3)
format(rpluy(y=recent$chpermwage,yo=schpermwage1,
             ywgt=recent$wgt,yowgt=wschpermwage1,
             upper=TRUE, pvalue=TRUE),
       digits=3)
#
# original$chpermwage is the change in perm wages for the Original
# cohort (with nobs > 2 and non-attrited)
#
# el1 is the final education for the Original
# cohort (with nobs > 2 and non-attrited)
# 1= < HS 2= HS 3=HS+ 4= College+
#
# First generate the samples
#
#
el1 <- e1
el1[el1 < 12] <- 1
el1[el1 == 12] <- 2
el1[el1 > 12 & el1 < 16] <- 3
el1[el1 >= 16] <- 4
#
el2 <- e2
el2[el2 < 12] <- 1
el2[el2 == 12] <- 2
el2[el2 > 12 & el2 < 16] <- 3
el2[el2 >= 16] <- 4
#
i3x <- sample(seq(along=original$chpermwage), 
            size = 10*length(original$chpermwage), 
            prob=rdsamp(el2,el1,recent$wgt,original$wgt),
		 replace = TRUE)
schpermwage1 <- original$chpermwage[i3x]
wschpermwage1 <- original$wgt[i3x]
sel1 <- el1[i3x]
#
pwhso  <- schpermwage1[sel1 <= 2 & !is.na(sel1)]
pwsco  <- schpermwage1[sel1 >  2 & !is.na(sel1)]
wgthso <- wschpermwage1[sel1 <= 2 & !is.na(sel1)]
wgtsco <- wschpermwage1[sel1 >  2 & !is.na(sel1)]
#
pwhsr  <- recent$chpermwage[el2 <= 2 & !is.na(el2)]
pwscr  <- recent$chpermwage[el2 >  2 & !is.na(el2)]
wgthsr <- recent$wgt[el2 <= 2 & !is.na(el2)]
wgtscr <- recent$wgt[el2 >  2 & !is.na(el2)]
#
postscript("Fig7.ps",width=5.5,height=5.5,horizontal=TRUE)
par(err=-1)
par(mfrow=c(2,2))
#
spwhso<-sample(pwhso,size=(100000),prob=wgthso/sum(wgthso),replace=TRUE)
spwsco<-sample(pwsco,size=(100000),prob=wgtsco/sum(wgtsco),replace=TRUE)
spwhsr<-sample(pwhsr,size=(100000),prob=wgthsr/sum(wgthsr),replace=TRUE)
spwscr<-sample(pwscr,size=(100000),prob=wgtscr/sum(wgtscr),replace=TRUE)
#
# Try to overlapping densities
#
nbar <- log(length(spwhso), base = 2) + 1
kwidth <- diff(range(spwhso))/nbar * 0.5
#kwidth <- 1.7*1.4*kwidth
kwidth <- 1.4*kwidth
#
dens1 <- density(spwhso, n = 500, width=2.0*kwidth)
plot(x = (dens1$x), y = dens1$y, type = "l", 
  xlab = "change in log permanent wage", ylab = "density", 
  xlim = c(-1, 3),
  ylim=c(0,1.2))
fig7legend <- list(x=c(0.9,0.9),y=c(1.25,1.25))                  
legend(fig7legend,lty=1:2,cex=0.5, bty="n",
legend=c("original cohort","recent cohort"))       
title(main=paste("(a) high-school or less"),cex=0.6)
#
dens2 <- density(spwhsr, n = 500, width=1.7*kwidth)
lines(x = (dens2$x), y = dens2$y, type = "l",lty=2) 
#
# First the two total r.d. of HS to C between cohort
#
g10hs <- reldist(y=pwhsr, yo=pwhso, ci=FALSE, smooth=0.4,
	ywgt=wgthsr, yowgt=wgthso,
	bar=TRUE, quiet=FALSE,
	ylim=c(0,4),
        xlab="proportion of the original cohort")
title(main=paste("(b) entropy = ",format(g10hs$entropy,digits=3)),cex=0.6)
abline(h=1,lty=2)
#
nbar <- log(length(spwscr), base = 2) + 1
kwidth <- diff(range(spwscr))/nbar * 0.5
#kwidth <- 1.75*1.2*kwidth
kwidth <- 1.2*kwidth
#
dens1 <- density(spwsco, n = 500, width=1.5*kwidth)
plot(x = (dens1$x), y = dens1$y, type = "l", 
  xlab = "change in log permanent wage", ylab = "density", 
  xlim = c(-1, 3),
  ylim=c(0,1.2))
fig7legend <- list(x=c(0.9,0.9),y=c(1.25,1.25))                  
legend(fig7legend,lty=1:2,cex=0.5, bty="n",
legend=c("original cohort","recent cohort"))       
title(main=paste("(c) more than high school"),cex=0.6)
#
dens2 <- density(spwscr, n = 500, width=2*kwidth)
lines(x = (dens2$x), y = dens2$y, type = "l",lty=2) 
#
g10sc <- reldist(y=pwscr, yo=pwsco, ci=FALSE, smooth=0.4,
	ywgt=wgtscr, yowgt=wgtsco,
	bar=TRUE, quiet=FALSE,
	ylim=c(0,4),
        xlab="proportion of the original cohort")
title(main=paste("(d) entropy = ",format(g10sc$entropy,digits=3)),cex=0.6)
abline(h=1,lty=2)
#
format(100*cbind(g10hs$deciles,g10sc$deciles,
                 g10sc$deciles-g10hs$deciles),digits=3)
#
collectresults <- c(
  wtd.quantile(pwhsr,weight=wgthsr),
  wtd.quantile(pwhso,weight=wgthso),
  wtd.quantile(pwhsr,weight=wgthsr)-wtd.quantile(pwhso,weight=wgthso),
  exp(wtd.quantile(pwhsr,weight=wgthsr)-wtd.quantile(pwhso,weight=wgthso)),
  g10hs$entropy,
  g10hs$entropy-reldist(y=pwhsr, yo=pwhso, ywgt=wgthsr, yowgt=wgthso, quiet=FALSE,
                     show="residual",graph=FALSE)$entropy,
  reldist(y=pwhsr, yo=pwhso, ywgt=wgthsr, yowgt=wgthso, quiet=FALSE,
                     show="residual",graph=FALSE)$entropy,
  g10hs$rp[2],
  g10hs$rpl[2],
  g10hs$rpu[2]
)
#
collectresults <- cbind(collectresults,c(
  wtd.quantile(pwscr,weight=wgtscr),
  wtd.quantile(pwsco,weight=wgtsco),
  wtd.quantile(pwscr,weight=wgtscr)-wtd.quantile(pwsco,weight=wgtsco),
  exp(wtd.quantile(pwscr,weight=wgtscr)-wtd.quantile(pwsco,weight=wgtsco)),
  g10sc$entropy,
  g10sc$entropy-reldist(y=pwscr, yo=pwsco, ywgt=wgtscr, yowgt=wgtsco, quiet=FALSE,
                     show="residual",graph=FALSE)$entropy,
  reldist(y=pwscr, yo=pwsco, ywgt=wgtscr, yowgt=wgtsco, quiet=FALSE,
                     show="residual",graph=FALSE)$entropy,
  g10sc$rp[2],
  g10sc$rpl[2],
  g10sc$rpu[2])
)
#
dimnames(collectresults)<-list(
  c("median log recent","median log original",
    "diff median log", "median ratio (recent/original)",
    "entropy", "median entropy","shape entropy","MRP","LRP","URP"),
  c("high-school or less","some college")
)
#
print(collectresults)
binn <- 10
postscript(file = "Fig8.ps",width=6.5,height=3.5,horizontal=TRUE)
par(err=-1)
par(mfrow=c(1,2))
#
rdhsrscr <- rdeciles(y=pwhsr, yo=pwscr, ywgt=wgthsr, yowgt=wgtscr, binn=binn)
rdhsosco <- rdeciles(y=pwhso, yo=pwsco, ywgt=wgthso, yowgt=wgtsco, binn=binn)
#
mscrdhsrscr <- rdeciles(
 y=pwhsr - wtd.quantile(pwhsr, weight=wgthsr) + wtd.quantile(pwhso, weight=wgthso),
yo=pwscr - wtd.quantile(pwscr, weight=wgtscr) + wtd.quantile(pwsco, weight=wgtsco),
ywgt=wgthsr, yowgt=wgtscr, binn=binn)
#
mhsrdhsrscr <- rdeciles(
 y=pwhso - wtd.quantile(pwhso, weight=wgthso) + wtd.quantile(pwhsr, weight=wgthsr),
yo=pwsco - wtd.quantile(pwsco, weight=wgtsco) + wtd.quantile(pwscr, weight=wgtscr),
ywgt=wgthso, yowgt=wgtsco, binn=binn)
#
m1rdhsrscr <- rdeciles(yo=pwsco, 
y=pwhsr - wtd.quantile(pwhsr, weight=wgthsr) + wtd.quantile(pwhso, weight=wgthso),
yowgt=wgtsco, ywgt=wgthsr, binn=binn)
#
m2rdhsrscr <- rdeciles(y=pwhso, 
yo=pwscr - wtd.quantile(pwscr, weight=wgtscr) + wtd.quantile(pwsco, weight=wgtsco),
ywgt=wgthso, yowgt=wgtscr, binn=binn)
#
m3rdhsrscr <- rdeciles(y=pwhsr, 
yo=pwsco - wtd.quantile(pwsco, weight=wgtsco) + wtd.quantile(pwscr, weight=wgtscr),
yowgt=wgtsco, ywgt=wgtscr, binn=binn)
#
achange <- binn*(rdhsrscr$x - rdhsosco$x)
armeff  <- binn*(mhsrdhsrscr$x - rdhsosco$x)
ahseff  <- binn*(m1rdhsrscr$x - rdhsosco$x)
asceff  <- binn*(m2rdhsrscr$x - rdhsosco$x)
ainteff <- achange - armeff - ahseff - asceff
#
bloc <- (1:binn)-0.5
barplot(height=achange,space=0,width=1,axes=FALSE,
     xlab="Decile",ylab="Percentage Point Change",
     ylim=c(-20.0,25))
axis(1,labels=paste(1:binn),at=bloc)
axis(2,labels=TRUE,at=seq(-20.0,25,length=10))
title(main="(a) Marginal effects",cex=0.6)
lines(y=(armeff),x=bloc,lty=1)
lines(y=(asceff),x=bloc,lty=3)
lines(y=(ahseff),x=bloc,lty=2)
abline(h=seq(-20,25,length=10),lty=2)
points(y=(armeff),x=bloc,pch=16,cex=0.7)
points(y=(asceff),x=bloc,pch=3,cex=0.7)
points(y=(ahseff),x=bloc,pch=1,cex=0.7)
#
fig8legend <- list(x=c(4,4),y=c(24.3,24.3))                  
legend(fig8legend,pch=c(16,1,3),lty=c(1:3),cex=0.5, bty="n",
legend=c("Change in relative median",
         "High-school shape effect","College shape effect"))
#
armeff  <- binn*(mhsrdhsrscr$x - rdhsosco$x)
ahseff  <- binn*(m3rdhsrscr$x - mhsrdhsrscr$x)
asceff  <- binn*(rdhsrscr$x - m3rdhsrscr$x)
#
barplot(height=achange,space=0,width=1,axes=FALSE,
     xlab="Decile",ylab="Percentage Point Change",
     ylim=c(-20.0,25))
axis(1,labels=paste(1:binn),at=bloc)
axis(2,labels=TRUE,at=seq(-20.0,25,length=10))
title(main="(b) Sequential effects",cex=0.6)
lines(y=(armeff),x=bloc,lty=1)
lines(y=(asceff),x=bloc,lty=3)
lines(y=(ahseff),x=bloc,lty=2)
abline(h=seq(-20,25,length=10),lty=2)
points(y=(armeff),x=bloc,pch=16,cex=0.7) # filled octagon
points(y=(asceff),x=bloc,pch=3,cex=0.7)  # cross
points(y=(ahseff),x=bloc,pch=1,cex=0.7)  # octagon
#
fig8legend <- list(x=c(4,4),y=c(24.3,24.3))                  
legend(fig8legend,pch=c(16,1,3),lty=c(1:3),cex=0.5, bty="n",
legend=c("Change in relative median",
         "High-school shape effect","College shape effect"))
