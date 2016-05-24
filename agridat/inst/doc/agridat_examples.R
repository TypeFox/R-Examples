## ----setup, echo=FALSE, results="hide"--------------------------------------------------
library("knitr")
opts_chunk$set(fig.align="center", fig.width=7, fig.height=7)
options(width=90)

## ----packs, message=FALSE---------------------------------------------------------------
library("agridat")
library("HH")
library("lattice")
library("latticeExtra")
library("mapproj")
library("maps")
library("reshape2")

## ----lee1, echo=FALSE, fig.height=7.5---------------------------------------------------
dat <- lee.potatoblight
# Note the progression to lower scores as time passes in each year
skp <- c(rep(0,10),
         rep(0,7),1,1,1,
         rep(0,8),1,1,
         rep(0,6),1,1,1,1,
         rep(0,5),1,1,1,1,1,
         rep(0,5),1,1,1,1,1,
         rep(0,6),1,1,1,1,
         rep(0,5),1,1,1,1,1,
         rep(0,5),1,1,1,1,1,
         rep(0,5),1,1,1,1,1)
desplot(y ~ col*row|date, dat,
        main="lee.potatoblight",
        between=list(y=.3), strip.cex =.7,
        layout=c(10,11), skip=as.logical(skp))

## ----lee2, echo=FALSE-------------------------------------------------------------------
# 1983 only.  I.Hardy succumbs quickly
dat <- lee.potatoblight
dat$dd <- as.Date(dat$date)
d83 <- droplevels(subset(dat, year==1983))
foo <- xyplot(y ~ dd|gen, d83, group=rep,
       xlab="Date", ylab="Blight resistance score",
       main="lee.potatoblight 1983", as.table=TRUE,
       par.settings=list(
         superpose.symbol=list(col=c("black","red","royalblue","#009900","dark orange"),
           pch=c("1","2","3","4","5"))),
       scales=list(alternating=FALSE, x=list(rot=90, cex=.7)))
foo + xyplot(y ~ dd|gen, d83, subset=year==1983, type='smooth', col='gray80')

## ----harrison, echo=FALSE, fig.height=6-------------------------------------------------
d1 <- subset(harrison.priors, substance=="daidzein")
d1 <- d1[ , c("source","number","min","max")]
out <- data.frame(source=vector("character"),
                  vals=vector("numeric"))
for(ii in 1:nrow(d1)){
  n <- d1[ii,'number']
  mi <- d1[ii,'min']; ma <- d1[ii,'max']
  mu <- mean(log(c(mi,ma)))
  sig <- (log(mi) - mu ) / qnorm(1/(1+n))
  vals <- exp(mu + sig*qnorm(seq(1/(1+n), to=n/(1+n), length=n)))
  out <- rbind(out, data.frame(source=d1[ii,'source'], vals=vals))
}
out <- droplevels(out) # Extra levels exist in d1
foo0 <- dotplot(source ~ vals, out, main="harrison.priors", xlab="Daidzein level",
                panel=function(x,y,...){
                  panel.dotplot(x,y,...)
                  #browser()
                  # Minimum for each row
                  x2l <- tapply(x, y, min)
                  x2r <- tapply(x, y, max)
                  y2 <- tapply(y, y, "[", 1)
                  panel.xyplot(x2l, y2, pch=16, cex=1.5, col="navy")
                  panel.xyplot(x2r, y2, pch=16, cex=1.5, col="navy")
                },
                # Hack.  Add blanks for extra space on graph
                ylim=c(levels(out$source),"","","","prior","Constructed","",""))

# Now calculate parameters for a common lognormal distribution
mu0 <- mean(log(out$vals))
sd0 <- sd(log(out$vals))
xvals <- seq(0,2000, length=100)
library(latticeExtra)
foo0 + xyplot((19+4000*dlnorm(xvals, mu0, sd0))~xvals, type='l',
              panel=function(x,y,...){
                panel.xyplot(x,y,...)
                panel.abline(h=19, col="gray90")
              })


## ----mead, echo=FALSE-------------------------------------------------------------------
dat <- mead.germination
# dat <- transform(dat, concf=factor(conc))

# Use log conc as a covariate.

# Note, my graph showing the density of the fit is similar to graphs
# found at the following site.  Did I get my idea here?  I don't know.
# http://www.unc.edu/courses/2010fall/ecol/563/001/docs/lectures/lecture18.htm#numerical

dat <- transform(dat, logconc=log(conc+.01))
m6 <- glm(cbind(germ, seeds-germ) ~ temp + temp:logconc, data=dat,
          family=binomial(link="logit"))

# Estimates of p for the binomial densities
newb <- expand.grid(temp=c('T1','T2','T3','T4'), logconc=log(c(0,.1,1,10)+.01))
newb$pct <- predict(m6, new=newb, type='response')
# Binomial density
foob <- xyplot(pct~logconc |temp, newb,
               xlim=c(-5.5, 4.5), ylim=c(-2, 53), as.table=TRUE,
               xlab="Log concentration",
               ylab="Seeds germinating (out of 50).  Binomial density.",
               main="mead.germination", #layout=c(4,1),
               panel=function(x,y,...){
                 for(ix in 1:4){
                   quan <- qbinom(c(.025, .975), size=50, prob=y[ix])
                   yval <- seq(min(quan), max(quan), by=1)
                   off <- x[ix]
                   xl <- off + rep(0, length(yval))
                   # Constant multiuplier of 8 chosen by trial and error
                   xr <- off + 8 * dbinom(yval, size=50, prob=y[ix])
                   panel.segments(xl,yval,xr, yval, cex=.35, lwd=3, col="gray70")
                 }
               })


# Add mean response line with equally-spaced points on the log scale
newl <- expand.grid(temp=c('T1','T2','T3','T4'),
                   logconc=seq(log(.01), log(10.01), length=50))
newl$pct <- predict(m6, new=newl, type='response')
# Logistic curve
fool <- xyplot(pct~logconc|temp, newl,
               panel=function(x,y,...){
                 panel.points(x, 50*y, type='l', col='blue')
               })


# Data points last, on top of everything
food <- xyplot(germ~logconc|temp, dat, layout=c(4,1),
       ylab="Seeds germinating (out of 50)", cex=1.5, pch=20, col='black')
foob + fool + food


## ----gomez, echo=FALSE------------------------------------------------------------------
dat <- gomez.stripsplitplot

# Layout
desplot(gen~x+y, dat, out1=rep, col=nitro, text=planting, cex=1,
        main="gomez.stripsplitplot")

## ----gomez2, echo=FALSE-----------------------------------------------------------------
dat <- gomez.splitsplit
dat$nitrogen <- factor(dat$nitro)
if(require(HH)){
  position(dat$nitrogen) <- c(0,50,80,110,140)
  interaction2wt(yield~rep+nitrogen+management+gen, data=dat,
                 main="gomez.splitsplit",
                 relation=list(x="free", y="same"),
                 rot=c(90,0), xlab="",
                 par.strip.text.input=list(cex=.8))
}

## ----keen, echo=FALSE, fig.width=7, fig.height=7.5--------------------------------------
dat <- keen.potatodamage

# Energy E1, Rod R4, Weight W1 have higher proportions of severe damage
# Rod 8 has the least damage
d2 <- xtabs(count~energy+rod+gen+weight+damage, data=dat)
mosaicplot(d2, color=c("lemonchiffon1","moccasin","lightsalmon1","indianred"),
           xlab="Energy / Genotype", ylab="Rod / Weight", main="keen.potatodamage",
           off=c(3,10,10,8,0),border="gray50")


## ----wright, echo=FALSE-----------------------------------------------------------------
dat <- minnesota.barley.yield
datw <- minnesota.barley.weather

# Weather trends over time
library(latticeExtra)
#useOuterStrips(xyplot(cdd~mo|year*site, datw, groups=year,
#main="minnesota.barley", xlab="month", ylab="Cooling degree days",
#subset=(mo > 3 & mo < 10), scales=list(alternating=FALSE),
#type='l', auto.key=list(columns=5)))

# Total cooling/heating/precip in Apr-Aug for each site/yr
ww <- subset(datw, mo>=4 & mo<=8)
ww <- aggregate(cbind(cdd,hdd,precip)~site+year, data=ww, sum)

# Average yield per each site/env
yy <- aggregate(yield~site+year, dat, mean)

minn <- merge(ww, yy)


# Higher yields generally associated with cooler temps, more precip
library(reshape2)
me <- melt(minn, id.var=c('site','year'))
mey <- subset(me, variable=="yield")
mey <- mey[,c('site','year','value')]
names(mey) <- c('site','year','y')
mec <- subset(me, variable!="yield")
names(mec) <- c('site','year','covar','x')
mecy <- merge(mec, mey)
mecy$yr <- factor(mecy$year)
oldpar <- tpg <- trellis.par.get()
tpg$superpose.symbol$pch <- substring(levels(mecy$yr),4) # Last digit of year
trellis.par.set(tpg)
foo <- xyplot(y~x|covar*site, data=mecy, groups=yr, cex=1, ylim=c(5,65),
              xlab="Weather covariate", ylab="Barley yield",
              main="minnesota.barley",
              panel=function(x,y,...) {
                panel.lmline(x,y,..., col="gray")
                panel.superpose(x,y,...)
              },
              scales=list(x=list(relation="free")))
foo <- useOuterStrips(foo, strip.left = strip.custom(par.strip.text=list(cex=.7)))
combineLimits(foo, margin.x=2L)

## ----crossa, echo=FALSE, message=FALSE--------------------------------------------------

# Specify env.group as column in data frame
dat2 <- crossa.wheat
dat2$eg <- ifelse(is.element(dat2$loc,
c("KN","NB","PA","BJ","IL","TC", "JM","PI","AS","ID","SC","SS",
"SJ","MS","MG","MM")), "Grp1", "Grp2")
m4 <- gge(yield~gen*loc, dat2, env.group=eg, scale=FALSE)
# plot(m4)
biplot(m4, lab.env=TRUE, title="crossa.wheat")

## ----nebr1, echo=FALSE------------------------------------------------------------------
library("maps")
library("mapproj")
library("latticeExtra")

dat <- nebraska.farmincome
dat$stco <- paste0('nebraska,', dat$county)
dat <- transform(dat, crop=crop/1000, animal=animal/1000)

# Raw, county-wide incomes.  Note the outlier Cuming county
mapplot(stco ~ crop + animal, data = dat,
        scales = list(draw = FALSE),
        main="nebraska.farmincome",
        xlab="", ylab="Income ($1000) per county",
        colramp=RedGrayBlue,
        map = map('county', 'nebraska', plot = FALSE, fill = TRUE,
                  projection = "mercator"))


## ----nebr2, echo=FALSE------------------------------------------------------------------

# Now scale to income/mile^2
dat <- within(dat, {
  crop.rate <- crop/area
  animal.rate <- animal/area
})
# And use manual breakpoints.
mapplot(stco ~ crop.rate + animal.rate, data = dat,
        scales = list(draw = FALSE),
        main="nebraska.farmincome",
        xlab="", ylab="Income ($1000) per square mile (percentile breaks)",
        map = map('county', 'nebraska', plot = FALSE, fill = TRUE,
                  projection = "mercator"),
        colramp=RedGrayBlue,
        #breaks=quantile(c(dat$crop.rate, dat$animal.rate),
        #                c(0,.1,.2,.4,.6,.8,.9,1), na.rm=TRUE)
        # To eliminate dependency on classInt package, hardcode the breakpoints
        #breaks=classIntervals(na.omit(c(dat$crop.rate, dat$animal.rate)), n=7, style='fisher')$brks
        breaks=c(0,.049, .108, .178, .230, .519, .958, 1.31)
        )

## ----lasrosas,echo=FALSE, fig.height=7.5------------------------------------------------

dat <- lasrosas.corn
library(latticeExtra)

# yield map
foo1 <- levelplot(yield ~ long*lat|factor(year), data=dat,
          aspect=1,
          main="lasrosas.corn grain yield (qu/ha)", xlab="Longitude", ylab="Latitude",
          scales=list(alternating=FALSE),
          prepanel = prepanel.default.xyplot,
          panel = panel.levelplot.points,
          type = c("p", "g"), col.regions=RedGrayBlue)

# Experiment design...shows problems in 2001
dat <- lasrosas.corn

xl <- range(dat$long)
yl <- range(dat$lat)

sseq=matrix(c(
  35, 0.9, 0.5,  # brown
  35, 0.8, 0.6,
  35, 0.7, 0.7,
  35, 0.6, 0.8,
  35, 0.5, .9,
  35, 0.4, 1,
  80, 0.9, 0.5,  # green
  80, 0.8, 0.6,
  80, 0.7, 0.7,
  80, 0.6, 0.8,
  80, 0.5, 0.9,
  80, 0.4, 1,
  190, 0.9, 0.5,  # blue
  190, 0.8, 0.6,
  190, 0.7, 0.7,
  190, 0.6, 0.8,
  190, 0.5, 0.9,
  190, 0.4, 1
  ), ncol=3, byrow=TRUE)
sseq <- hsv(sseq[,1]/360, sseq[,2], sseq[,3])

dat$repnf <- factor(paste(dat$rep,dat$nf))
# levels(dat$repnf) # check the order
#dat <- transform(dat, col=as.character(sseq[as.numeric(factor(paste(dat$rep,dat$nf)))]))

# By default, manual specification of col/pch does not work with multiple panels.
# Define a custom panel function to make it work
mypanel <- function(x,y,...,subscripts,col,pch) {
  panel.xyplot(x,y,col=col[subscripts],pch=pch[subscripts], ...)
}

foo2 <- xyplot(lat~long|factor(year), data=dat,
       aspect=1, xlim=xl, ylim=yl, cex=0.9,
       main="lasrosas.corn experiment design", xlab="", ylab="",
       scales=list(alternating=FALSE),
       col=sseq[dat$repnf],
       #pch=levels(dat$topo)[dat$topo],
       pch=c('-','+','/','\\')[dat$topo],
       panel=mypanel)

plot(foo1, split = c(1, 1, 1, 2))
plot(foo2, split = c(1, 2, 1, 2), newpage = FALSE)


## ----nass, echo=FALSE, fig.height=8-----------------------------------------------------
dat <- nass.corn
dat$acres <- dat$acres/1000000

# Use only states that grew at least 100K acres of corn in 2011
keep <- droplevels(subset(dat, year == 2011 & acres > .1))$state
dat <- subset(dat, state != "Delaware")
dat <- subset(dat, state != "Idaho")
dat <- subset(dat, state != "Washington")
dat <- subset(dat, state != "California")
dat <- droplevels(subset(dat, is.element(state, keep)))
# Acres of corn grown each year
xyplot(acres ~ year|state, dat, type='l', as.table=TRUE,
       layout=c(6,5),
       strip=strip.custom(par.strip.text=list(cex=.5)),
       main="nass.corn", xlab="Year", ylab="Million acres of corn")


## ----finish, echo=FALSE, results="asis"-------------------------------------------------
toLatex(sessionInfo(), locale=FALSE)

