# course.r
#
library(gamlss)
library(AGD)
boys <- boys7482

setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
old.par <- par(mar=c(1,1,1,1), col="darkblue")
dir  <- path.expand("~/Documents/Sync/Groeistat/NIHES/2012/Slides/Session1Figures/")


# scatterplot hgt - age
postscript(paste(dir,"hgtage1.eps",sep=""))
par(mar=c(2,2,2,2), col="darkblue")
with(boys, plot(age, hgt))
dev.off()

# simple linear model
fit <- lm(hgt~age, data=boys)
summary(fit)

# plot diagnostics
postscript(paste(dir,"fit1_%01d.eps",sep=""))
par(mar=c(2,2,2,2), col="darkblue")
plot(fit)
dev.off()

# add quadratic
fit <- lm(hgt~age+I(age^2), data=boys, na.action=na.exclude)
summary(fit)

postscript(paste(dir,"fit2_%01d.eps",sep=""))
par(mar=c(2,2,2,2), col="darkblue")
plot(fit)
dev.off()

postscript(paste(dir,"hgtage2.eps",sep=""))
par(mar=c(2,2,2,2), col="darkblue")
with(boys, plot(age, hgt))
sigma <- summary(fit)$sigma
lines(boys$age, predict(fit), col="red")
lines(boys$age, predict(fit)+2*sigma, col="green")
lines(boys$age, predict(fit)-2*sigma, col="green")
dev.off()

# semi-parametric model
inner <- c(7/365.25, 1/3, 1, 2, 4, 10, 14, 18, 21)
xb <- bs(boys$age, knots=inner, B=c(0,23), deg=1)
xb[1:3,]

fit <- lm(hgt ~ xb, data=boys, na.action=na.exclude)

postscript(paste(dir,"bspline1.eps",sep=""))
par(mar=c(4,4,1,1), col="darkblue")
with(boys, plot(age, hgt))
lines(boys$age, predict(fit), col="red", lwd=3)
dev.off()

# include fitted knots
postscript(paste(dir,"bspline2.eps",sep=""))
par(mar=c(4,4,1,1), col="darkblue")
with(boys, plot(age, hgt))
lines(boys$age, predict(fit), col="red", lwd=3)
est <- coef(fit)
xy <- xy.coords(x=inner[1:9], y=(est[1]+est)[c(-1,-11)])
points(xy, col="yellow", pch=20, cex=2.5)
dev.off()

# include +2SD and -2SD lines
postscript(paste(dir,"bspline3.eps",sep=""))
par(mar=c(4,4,1,1), col="darkblue")
with(boys, plot(age, hgt))
lines(boys$age, predict(fit), col="red", lwd=3)
est <- coef(fit)
xy <- xy.coords(x=inner[1:9], y=(est[1]+est)[c(-1,-11)])
points(xy, col="yellow", pch=20, cex=2.5)
sigma <- summary(fit)$sigma
lines(boys$age, predict(fit)+2*sigma, col="green", lwd=3)
lines(boys$age, predict(fit)-2*sigma, col="green", lwd=3)
dev.off()

myref <- data.frame(sex   = "M",
                    sub   = "N",
                    x   = round(xy$x,4),
                    mean  = round(xy$y,2),
                    sd    = round(sigma,2))

# calculate SDS of length 80cm, 1 years old boy
# relative to your new reference
y2z(y=80, x=1, ref=myref, dist="NO")

# relative to Dutch reference (default)
y2z(y=80,x=1)

# relative to CDC reference
y2z(y=80,x=1,ref=cdc.hgt)

# calculate PERCENTILE of length 80cm, 1 years old boy
# relative to your new reference
100*pnorm(y2z(y=80,x=1,ref=myref,dist="NO"))

# relative to Dutch reference
100*pnorm(y2z(y=80,x=1))

# relative to CDC reference
100*pnorm(y2z(y=80,x=1,ref=cdc.hgt))

# calculate PERCENTILE of length 80cm, 1 years old GIRL
# relative to your new reference
# The warning indicates that there are no female references
100*pnorm(y2z(y=80,x=1,sex="F",ref=myref,dist="NO"))

# relative to Dutch reference
100*pnorm(y2z(y=80,x=1,sex="F",ref=nl4.hgt))

# relative to CDC reference
100*pnorm(y2z(y=80,x=1,sex="F",ref=cdc.hgt))

# SDS of IOTF BMI cut-off value for overweight (boys 2-18)
# relative to Dutch boys reference
cutoff <- c(
18.41, 18.15, 17.89, 17.72, 17.55, 17.49, 17.42, 17.49, 17.55, 17.74,
17.92, 18.18, 18.44, 18.77, 19.10, 19.47, 19.84, 20.20, 20.55, 20.89,
21.22, 21.57, 21.91, 22.27, 22.62, 22.96, 23.29, 23.60, 23.90, 24.18,
24.46, 24.73, 25.00)
age <- seq(2, 18, by=0.5)
percent <- 100-100*pnorm(y2z(y=cutoff, x=age, sex="M", ref=nl4.bmi))

postscript(paste(dir,"iotf1.eps",sep=""))
par(mar=c(4,4,1,1), col="darkblue")
plot(age, percent, type='l',ylim=c(0,20),lwd=2,col="red")
title("Overweight prevalence Dutch boys 1997")
dev.off()

# percentage of children lighter than 15kg at ages 2-5
e <- expand.grid(age=2:5, sex=c("M","F"))
z <- y2z(y=rep(15,nrow(e)), x=e$age, sex=e$sex, ref=nl4.wgt)
w <- matrix(100*round(pnorm(z),2), nrow=2, byrow=TRUE)
dimnames(w) <- list(c("boys","girls"),2:5)
w

# analysis in Z scale
hgt.z <- y2z(boys$hgt, boys$age, sex="M", ref=nl4.hgt)
wgt.z <- y2z(boys$wgt, boys$age, sex="M", ref=nl4.wgt)
postscript(paste(dir,"hgtwgt1.eps",sep=""))
par(mar=c(4,4,1,1), col="darkblue")
plot(hgt.z, wgt.z)
dev.off()

# standard set of Z-scores of weight for all tabulated ages, boys & girls
sds <- c(-2.5, -2, -1, 0, 1, 2, 2.5)
age <- nl4.wgt$x
z <- rep(sds, times=length(age))
x <- rep(age, each=length(sds))
sex <- rep(c("M","F"), each=length(z)/2)
w <- z2y(z=z, x=x, sex=sex, ref=nl4.wgt)
w <- matrix(w, ncol=length(sds), byrow=TRUE)
dimnames(w) <- list(age, sds)
postscript(paste(dir,"wfa1.eps",sep=""))
par(mar=c(4,4,1,1), col="darkblue")
matplot(x=age[1:69], y=w[1:69,], xlim=c(0,21),
       type="l", lty=1, lwd=2,
       col=c("red","green","gray","blue","gray","green","red"),
       xlab="Age (years)", ylab="Weight (kg)")
dev.off()

# interpolate standard to days
days <- 0:15
sds  <- c(-2, 0, +2)
z    <- rep(sds, length(days))
x    <- rep(round(days/365.25,4), each=length(sds))
w    <- z2y(z, x, sex="M", ref=nl4.hgt)
w    <- matrix(w, ncol=length(sds), byrow=TRUE)
dimnames(w) <- list(days, sds)
w



# --- SESSION 3
# --- DISTRIBUTIONS  

par(old.par)
par(cex=1.5, lwd=1.2, col="darkblue")
dir  <- path.expand("~/Documents/Sync/Groeistat/NIHES/2012/Slides/Session3Figures/")
lwd <- 1.5
cex <- 1.2

# boys distribution 10 years olds
nl4.hgt[46:48,]
m  <- nl4.hgt[47,"M"]
m
sd <- nl4.hgt[47,"S"]*m
sd

# normal distribution
y <- 120:170
d <- dNO(y, mu=m, sigma=sd)   # NOTE: use dNO, alternative to dnorm
postscript(paste(dir,"f_dnorm.eps",sep=""))
par(lwd=lwd, cex=cex)
plot(y, d, type = 'l', col = 'blue', 
     xlab="Measurement",ylab="Density")
dev.off()

# cumulative normal distribution
y <- 120:170
p <- pNO(y, mu=m, sigma=sd)   
postscript(paste(dir,"f_pnorm.eps",sep=""))
par(lwd=lwd, cex=cex)
plot(y, p, type = 'l', col = 'blue', 
     xlab="Measurement", ylab="Probability")
dev.off()

# normal quantiles
p <- seq(0.001, 0.999, by =0.001)
y <- qNO(p)   # NOTE: uses qNO of gamlss instead of qnorm
postscript(paste(dir,"f_qnorm.eps",sep=""))
par(lwd=lwd, cex=cex)
plot(p, y, type = 'l', col = 'blue', 
      xlab="Probability", ylab="Measurement")
dev.off()

# with data
boys <- boys7482
y <- boys[boys$age>9.5 & boys$age<10.5,"hgt"]
postscript(paste(dir,"f_hnorm1.eps",sep=""))
par(lwd=lwd, cex=cex)
hist(y, xlab="Height (cm)", main="")
dev.off()

# play with histogram
postscript(paste(dir,"f_hnorm2.eps",sep=""))
par(lwd=lwd, cex=cex)
b <- seq(120,166,2)
hist(y, xlab="Height (cm)", main="", breaks=b, col='gray', border='white')
dev.off()

# add normal references
postscript(paste(dir,"f_hnorm3.eps",sep=""))
par(lwd=lwd, cex=cex)
hist(y, xlab="Height (cm)", main="", breaks=b, col='gray', border='white',freq=F)
lines(120:166, dNO(120:166,mu=m,sig=sd), col = 'blue', lwd = 2)
dev.off()

# 
postscript(paste(dir,"f_hnorm4.eps",sep=""))
par(lwd=lwd, cex=cex)
b <- seq(120,166,2)
hist(y, xlab="Height (cm)", main="", breaks=b, col='gray', border='white',freq=F)
lines(120:166, dNO(120:166,mu=m,sig=sd), col = 'blue')
lines(120:166, dNO(120:166,mu=140.8,sig=0.0454*140.8), col = 'red', lty=2)
lines(120:166, dNO(120:166,mu=145.7,sig=0.0465*145.7), col = 'red', lty=2)
dev.off()

# empirical cumulative distribution
y <- boys[boys$age>9.5 & boys$age<10.5,"hgt"]
n <- length(y)
p <- ((1:n) - 0.5) / n   # Empirical probabilities
postscript(paste(dir,"f_cumnorm.eps",sep=""))
par(lwd=lwd, cex=cex)
plot(sort(y), p, type = 's', col = 'red', xlab="Height (cm)")
dev.off()

# adding the theoretical cumulative distribution
n <- length(y)
p <- ((1:n) - 0.5) / n   # Empirical probabilities
postscript(paste(dir,"f_cumnorm2.eps",sep=""))
par(lwd=lwd, cex=cex)
plot(sort(y), p, type = 's', col = 'red', xlab="Height (cm)")
lines(120:170, pNO(120:170,m=m,sigma=sd), col = 'blue')
dev.off()


# Box-Cox transformation
y <- 20:70
postscript(paste(dir,"f_bc1.eps",sep=""))
par(lwd=lwd, cex=cex)
plot(function(y) dBCCG(y, mu=34,sigma=.16,nu=-1), 20, 70, ylab="Density", xlab="Weight (kg)", col="red")
dev.off()

# vary the nu parameter
postscript(paste(dir,"f_bc2.eps",sep=""))
par(lwd=lwd, cex=cex)
plot(function(y) dBCCG(y, mu=34,sigma=.16,nu=-1), 20, 70, ylab="Density", xlab="Weight (kg)", col="red")
lines(x=y, y=dBCCG(y, mu=34,sigma=.16,nu=0),col="green")
lines(x=y, y=dBCCG(y, mu=34,sigma=.16,nu=1),col="black")
lines(x=y, y=dBCCG(y, mu=34,sigma=.16,nu=2),col="blue")
dev.off()

# now with data
y <- boys[boys$age>9.5 & boys$age<10.5,"wgt"]
b <- seq(20,70,2)
postscript(paste(dir,"f_bc3.eps",sep=""))
par(lwd=lwd, cex=cex)
hist(y, xlab="Weight (kg)", main="", breaks=b, col='gray', border='white')
dev.off()

y <- sort(y)
# add 
postscript(paste(dir,"f_bc4.eps",sep=""))
par(lwd=lwd, cex=cex)
hist(y, xlab="Weight (kg)", main="", breaks=b, col='gray', border='white', freq=FALSE)
lines(x=y, y=dBCCG(y, mu=33.8,sigma=.162,nu=0.162),col="blue",lwd=2)
dev.off()





# LMS method
# BCCG
# BCT
# BCPE
# free distribution

# --- modelling four moments
# GAMLSS framework
# crosstable linear/smooth, homo/hetero
# linear model, homoskedastic (parametric)
# linear model, heteroskedastic (parametric)
# smooth model, homoskedastic
# smooth model, smooth heteroskedastic
# third and fourth moments

# --- SESSION 6
# --- DIAGNOSTICS
par(old.par)
par(cex=1.5, lwd=2, col="darkblue")
dir  <- path.expand("~/Documents/Sync/Groeistat/NIHES/2012/Slides/Session6Figures/")
lwd <- 2
cex <- 1.4

# if necessary first download file, and unzip to get to figures.R and the data
#err <- download.file(url="http://www.stefvanbuuren.nl/wormplot/code%20and%20data.zip", 
#        destfile="S:\\projecten\\a-i\\groeistat\\NIHES\\2010\\R\\download.zip")

library(gamlss)

# fitted.plot
data(abdom)
abd9 <- gamlss( y~cs(x,df=3), sigma.formula=~cs(x,df=3), 
                nu.formula=~1, tau.fomula=~1, family=BCT, data=abdom)
abd10 <- gamlss( y~cs(x,df=1), sigma.formula=~cs(x,df=1),
                nu.formula=~1, tau.fomula=~1, family=BCT, data=abdom)
postscript(paste(dir,"fitted.eps",sep=""))
par(lwd=1.7, cex.lab=2)
fittedPlot(abd9,abd10,x=abdom$x)
dev.off()

# plot points and centiles
postscript(paste(dir,"cent.eps",sep=""))
par(lwd=1.7, cex.lab=2)
centiles(abd9,xvar=abdom$x,xlab="Gestational age (week)",ylab="Circumference (cm)")
dev.off()

# empirical quantiles and fitted quantiles
# breaks <- nl4.hfa$age[1:69]

# calculate solution to demo worm plot
library(AGD)
data <- boys7482
data <- na.omit(data[,c("age","hgt","wgt")])
f0051  <- gamlss(hgt~cs(age,df=5,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(age,df=1,c.spar=c(-1.5,2.5)),
                 data=data,family=NO,
                 control=gamlss.control(n.cyc=3))
f0101  <- gamlss(hgt~cs(age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(age,df=1,c.spar=c(-1.5,2.5)),
                 data=data,family=NO,
                 control=gamlss.control(n.cyc=3))
f0151  <- gamlss(hgt~cs(age,df=15,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(age,df=1,c.spar=c(-1.5,2.5)),
                 data=data,family=NO,
                 control=gamlss.control(n.cyc=3))                      
f0154  <- gamlss(hgt~cs(age,df=15,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(age,df=4,c.spar=c(-1.5,2.5)),
                 data=data,family=NO,
                 control=gamlss.control(n.cyc=3))                      

# use transformed age
t.age <- fitted(lm(data$age~fitted(f0154)))
t.age <- t.age - min(t.age)
data.t <- data.frame(data,t.age=t.age)
f0051r <- gamlss(hgt~cs(t.age,df=5,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=1,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=2))
f0091r <- gamlss(hgt~cs(t.age,df=9,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=1,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=3))
f0101r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=1,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=2))
f0111r <- gamlss(hgt~cs(t.age,df=11,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=1,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=3))
f0105r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=5,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=3))
f0106r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=6,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=3))
f0107r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=7,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=3))
f1106r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=6,c.spar=c(-1.5,2.5)),
                 nu.formula=~cs(t.age,df=1,c.spar=c(-1.5,2.5)),
                 data=data.t,family=BCCG,
                 control=gamlss.control(n.cyc=3))
f4106r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=6,c.spar=c(-1.5,2.5)),
                 nu.formula=~cs(t.age,df=4,c.spar=c(-1.5,2.5)),
                 data=data.t,family=BCCG,
                 control=gamlss.control(n.cyc=3))



# utility functions
store <- function(..., path=.store){
    # stores object(s) in storage directory path (SvB May 2008)
    names <- as.character(substitute(list(...)))[-1]
    for (i in 1:length(names)){
        name <- names[i]
        file <- paste(path, name, sep="")
        file.access
        save(list=name,file=file)
    }
}                  

fetch <- function(..., path=.store){
    # fetches the object from the storage directory in path (SvB May 2008)
    names <- as.character(substitute(list(...)))[-1]
    for (i in 1:length(names)){
        name <- names[i]
        file <- paste(path, name, sep="")
        if (file.access(file,0)== -1) cat(paste("Warning: No file",file,"\n")) else load(file,.GlobalEnv)
    }
}

.store <- path.expand("~/Documents/Sync/Groeistat/NIHES/2010/R/store/")

# store(f0051,f0101,f0151,f0154,f0051r, f0091r, f0101r, 
#       f0111r, f0105r, f0106r, f0107r, f1106r, f4106r)
fetch(f0051,f0101,f0151,f0154,f0051r, f0091r, f0101r, 
      f0111r, f0105r, f0106r, f0107r, f1106r, f4106r)


postscript(paste(dir,"c0051.eps",sep=""))
centiles(f0051,xvar=data$age,main="0051",pch='.')
dev.off()
postscript(paste(dir,"c0101.eps",sep=""))
centiles(f0101,xvar=data$age,main="0101",pch='.')
dev.off()
postscript(paste(dir,"c0151.eps",sep=""))
centiles(f0151,xvar=data$age,main="0151",pch='.')
dev.off()
postscript(paste(dir,"c0051r.eps",sep=""))
centiles(f0051r,xvar=data$age,main="0051r",pch='.')
dev.off()
postscript(paste(dir,"c0106r.eps",sep=""))
centiles(f0106r,xvar=data$age,main="0106r",pch='.')
dev.off()
postscript(paste(dir,"c4106r.eps",sep=""))
centiles(f4106r,xvar=data$age,main="4106r",pch='.')
dev.off()

postscript(paste(dir,"w0051.eps",sep=""))
wp.twin(f0051,cex=0.5)
dev.off()
postscript(paste(dir,"w0101.eps",sep=""))
wp.twin(f0051,f0101,cex=0.5)
dev.off()
postscript(paste(dir,"w0151.eps",sep=""))
wp.twin(f0101,f0151,cex=0.5)
dev.off()
postscript(paste(dir,"w0051r.eps",sep=""))
wp.twin(f0051,f0051r,cex=0.5)
dev.off()
postscript(paste(dir,"w0101r.eps",sep=""))
wp.twin(f0051r,f0101r,cex=0.5)
dev.off()
postscript(paste(dir,"w0101r2.eps",sep=""))
wp.twin(f0091r,f0101r,cex=0.5)
dev.off()
postscript(paste(dir,"w0111r.eps",sep=""))
wp.twin(f0101r,f0111r,cex=0.5)
dev.off()
postscript(paste(dir,"w0106r.eps",sep=""))
wp.twin(f0101r,f0106r,cex=0.5)
dev.off()
postscript(paste(dir,"w0106r2.eps",sep=""))
wp.twin(f0105r,f0106r,cex=0.5)
dev.off()
postscript(paste(dir,"w0107r.eps",sep=""))
wp.twin(f0106r,f0107r,cex=0.5)
dev.off()
postscript(paste(dir,"w4106r.eps",sep=""))
wp.twin(f0106r,f4106r,cex=0.5)
dev.off()


# store(f0051,f0101,f0151,f0154,f0051r, f0091r, f0101r, 
#       f0111r, f0105r, f0106r, f0107r, f1106r, f4106r)
# 


#f0164 <- gamlss(hgt~cs(age,df=15), 
#               sigma.formula=~cs(age,df=4), 
#               data=data, family=NO, 
#               control=gamlss.control(n.cyc=3))

#par(lwd=1)
#wp.twin(fit3,cex=0.5)



# --- SESSION 8
# further tricks, outlook
par(cex=1.5, lwd=2, col="darkblue")
dir  <- path.expand("~/Documents/Sync/Groeistat/NIHES/2012/Slides/Session8Figures/")
lwd <- 2
cex <- 1.4

# age transformation
library(AGD)
data <- boys7482
data <- na.omit(data[,c("age","hgt","wgt")])
fit1 <- gamlss(hgt~cs(age,df=16), sigma.formula=~cs(age,df=4), data=data, family=NO, control=gamlss.control(n.cyc=5))
t.age <- fitted(lm(data$age~fitted(fit1)))
postscript(paste(dir,"agetra.eps",sep=""))
par(lwd=1.7, cex.lab=2)
plot(data$age, t.age, type='l',lwd=2,xlab="Age", ylab="Transformed age", col="blue",ylim=c(0,20),xlim=c(0,20))
abline(0,1,lwd=1)
dev.off()


# exporting a GAMLSS object to an LMS table
data2 <- cbind(data, t.age)
lms <- extractLMS(f0106r, data2)
tail(lms)

