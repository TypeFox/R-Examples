# R code for chapter 5 of Wood (2006) "GAMs: An Introduction with R"

## 5.1 Cherry trees again

library(mgcv)
data(trees)
ct1<-gam(Volume~s(Height)+s(Girth),
         family=Gamma(link=log),data=trees)
ct1
par(mfrow=c(1,2))
plot(ct1,residuals=TRUE)

## 5.1.1 Finer control of gam

ct2 <- gam(Volume~s(Height,bs="cr")+s(Girth,bs="cr"),
           family=Gamma(link=log),data=trees)
ct2

ct3 <- gam(Volume ~ s(Height)+s(Girth,bs="cr",k=20),
           family=Gamma(link=log),data=trees)
ct3

ct4 <- gam(Volume ~ s(Height) + s(Girth),
           family=Gamma(link=log),data=trees,gamma=1.4)
ct4
plot(ct4,residuals=TRUE)

## 5.1.2 Smooths of several variables

ct5 <- gam(Volume ~ s(Height,Girth,k=25),
           family=Gamma(link=log),data=trees)
ct5
plot(ct5,too.far=0.15)

ct6 <- gam(Volume ~ te(Height,Girth,k=5),
           family=Gamma(link=log),data=trees)
ct6
plot(ct6,too.far=0.15)

## 5.1.3 Parametric model terms

gam(Volume~Height+s(Girth),family=Gamma(link=log),data=trees)
trees$Hclass <- factor(floor(trees$Height/10)-5,
                labels=c("small","medium","large"))

ct7 <- gam(Volume ~ Hclass+s(Girth),
           family=Gamma(link=log),data=trees)
par(mfrow=c(1,2))
plot(ct7,all.terms=T)

anova(ct7)
AIC(ct7)
summary(ct7)

## 5.2 Brain imaging example
## 5.2.1 Preliminary modelling
library(gamair)
data(brain)
brain <- brain[brain$medFPQ>5e-3,] # exclude 2 outliers
m0 <- gam(medFPQ~s(Y,X,k=100),data=brain)
gam.check(m0) 

e <- residuals(m0); fv <- fitted(m0)
lm(log(e^2)~log(fv))

m1<-gam(medFPQ^.25~s(Y,X,k=100),data=brain)
gam.check(m1)

m2<-gam(medFPQ~s(Y,X,k=100),data=brain,family=Gamma(link=log),optimizer="perf")

mean(fitted(m1)^4);mean(fitted(m2));mean(brain$medFPQ)
m2

vis.gam(m2,plot.type="contour",too.far=0.03,
        color="gray",n.grid=60,zlim=c(-1,2))

## 5.2.2 Would an additive structure be better?
m3 <- gam(medFPQ~s(Y,k=30)+s(X,k=30),data=brain,
          family=Gamma(link=log),optimizer="perf")

m3
anova(m3,m2,test="F")

m4 <- gam(medFPQ~s(Y,k=30)+s(X,k=30)+s(Y,X,k=100),data=brain,
      family=Gamma(link=log),optimizer="perf")

## 5.2.3 Isotropic or tensor product smooths?

tm<-gam(medFPQ~te(Y,X,k=10),data=brain,family=Gamma(link=log),
        optimizer="perf")
tm1<-gam(medFPQ~s(Y,k=10,bs="cr")+s(X,bs="cr",k=10),
         data=brain,family=Gamma(link=log),optimizer="perf")

tm1
tm
anova(tm1,tm,test="F")

## 5.2.4 Detecting symmetry (with `by' variables)
brain$Xc <- abs(brain$X - 64.5)
brain$right <- as.numeric(brain$X<64.5)
m.sy <- gam(medFPQ~s(Y,Xc,k=100),data=brain,
            family=Gamma(link=log),optimizer="perf")
m.as <- gam(medFPQ~s(Y,Xc,k=100)+s(Y,Xc,k=100,by=right),
            data=brain,family=Gamma(link=log),optimizer="perf")
m.sy
m.as

anova(m.sy,m.as,test="F")
anova(m.as)

par(mfrow=c(1,3))
vis.gam(m.sy,plot.type="contour",view=c("Xc","Y"),too.far=.03,
        color="gray",n.grid=60,zlim=c(-1,2),main="both sides")
vis.gam(m.as,plot.type="contour",view=c("Xc","Y"),
        cond=list(right=0),too.far=.03,color="gray",n.grid=60,
        zlim=c(-1,2),main="left side")
vis.gam(m.as,plot.type="contour",view=c("Xc","Y"),
        cond=list(right=1),too.far=.03,color="gray",n.grid=60,
        zlim=c(-1,2),main="right side")

## 5.2.5 Comparing two surfaces

brain1 <- brain
mu <- fitted(m2)
n<-length(mu)
ind <- brain1$X<60 & brain1$Y<20
mu[ind] <- mu[ind]/3
set.seed(1)
brain1$medFPQ <- rgamma(rep(1,n),mu/m2$sig2,scale=m2$sig2)

brain2=rbind(brain,brain1)
brain2$sample1 <- c(rep(1,n),rep(0,n))
brain2$sample0 <- 1 - brain2$sample1

m.same<-gam(medFPQ~s(Y,X,k=100),data=brain2,
            family=Gamma(link=log),optimizer="perf")
m.diff<-gam(medFPQ~s(Y,X,k=100)+s(Y,X,by=sample1,k=100),
            data=brain2,family=Gamma(link=log),optimizer="perf")

anova(m.same,m.diff,test="F")

## 5.2.6 Prediction with `predict.gam'

predict(m2)[1:5]
pv <- predict(m2,se=TRUE)
pv$fit[1:5]
pv$se[1:5]

predict(m2,type="response")[1:5]
pv <- predict(m2,type="response",se=TRUE)
pv$se[1:5]

pd <- data.frame(X=c(80.1,68.3),Y=c(41.8,41.8))
predict(m2,newdata=pd)
predict(m2,newdata=pd,type="response",se=TRUE)

predict(m3,newdata=pd,type="terms",se=TRUE)

#### Prediction with `lpmatrix'

Xp <- predict(m2,newdata=pd,type="lpmatrix")
fv <- Xp%*%coef(m2)
fv

d <- t(c(1,-1))
d%*%fv
d%*%Xp%*%m2$Vp%*%t(Xp)%*%t(d)

## 5.2.7 Variances of non-linear functions of the fitted model

ind <- brain$region==1& ! is.na(brain$region)
Xp <- predict(m2,newdata=brain[ind,],type="lpmatrix")

library(MASS)
br <- mvrnorm(n=1000,coef(m2),m2$Vp) # simulate from posterior

mean.FPQ<-rep(0,1000)
for (i in 1:1000)
{ lp <- Xp%*%br[i,]  # replicate linear predictor
  mean.FPQ[i] <- mean(exp(lp)) # replicate region 1 mean FPQ
}

mean.FPQ <- colMeans(exp(Xp%*%t(br)))

hist(mean.FPQ)

## 5.3 Air pollution in Chicago example

data(chicago)
ap0 <- gam(death~s(time,bs="cr",k=200)+pm10median+so2median+
           o3median+tmpd,data=chicago,family=poisson)
gam.check(ap0)

par(mfrow=c(2,1))
plot(ap0,n=1000)  # n increased to make plot smooth
plot(ap0,residuals=TRUE,n=1000)

chicago$death[3111:3125]

ap1<-gam(death~s(time,bs="cr",k=200)+s(pm10median,bs="cr")+
     s(so2median,bs="cr")+s(o3median,bs="cr")+s(tmpd,bs="cr"),
     data=chicago,family=poisson)

lag.sum <- function(a,l0,l1)
## l0 is the smallest lag, l1 the largest
{ n<-length(a)
  b<-rep(0,n-l1)
  for (i in 0:(l1-l0)) b <- b + a[(i+1):(n-l1+i)]
  b
}
death<-chicago$death[4:5114]
time<-chicago$time[4:5114]
o3 <- lag.sum(chicago$o3median,0,3)
tmp <- lag.sum(chicago$tmpd,0,3)
pm10 <- lag.sum(log(chicago$pm10median+40),0,3)
so2 <- lag.sum(log(chicago$so2median+10),0,3)

ap2 <- gam(death ~ s(time,bs="cr",k=200) +
                   te(o3,tmp,pm10,k=c(8,8,6)),family=poisson)

ap3 <- gam(death ~ s(time,bs="cr",k=200) + te(o3,tmp,k=8) +
                   s(pm10,bs="cr",k=6),family=poisson)


## 5.4 Mackerel egg survey example
## 5.4.1 Model development

data(mack)
mack$log.net.area <- log(mack$net.area)

gm <- gam(egg.count~s(lon,lat,bs="ts")+s(I(b.depth^.5),bs="ts")
          +s(c.dist,bs="ts")+s(salinity,bs="ts")+s(temp.surf,bs="ts")
          +s(temp.20m,bs="ts")+offset(log.net.area),
          data=mack,family=poisson,scale=-1,gamma=1.4)

gm1 <- gam(egg.count ~ s(lon,lat,bs="ts",k=100)+
           s(I(b.depth^.5),bs="ts") + s(c.dist,bs="ts") +
           s(temp.surf,bs="ts") + s(temp.20m,bs="ts") +
           offset(log.net.area),
           data=mack,family=poisson,scale=-1,gamma=1.4)
gm1

gm1a <- gam(egg.count~s(lon,lat,bs="ts",k=100)+
            s(I(b.depth^.5),bs="ts") + s(temp.20m,bs="ts")+
            offset(log.net.area),data=mack,family=poisson,
            scale=-1,gamma=1.4)

gm2<-gam(egg.count ~ s(lon,lat,bs="ts",k=40) +
         s(I(b.depth^.5),bs="ts") + s(c.dist,bs="ts") +
         s(temp.surf,bs="ts") + s(temp.20m,bs="ts") +
         offset(log.net.area),
         data=mack,family=negbin(c(.01,100)),
         control=gam.control(maxit=100),gamma=1.4,optimizer="perf")
gm2

## alternative, using more recent methods...

gm2a<-gam(egg.count ~ s(lon,lat,bs="ts",k=40) +
         s(I(b.depth^.5),bs="ts") + s(c.dist,bs="ts") +
         s(temp.surf,bs="ts") + s(temp.20m,bs="ts") +
         offset(log.net.area),data=mack,family=negbin(c(.1,10)),
         control=gam.control(maxit=100),gamma=1.4)
gm2a

anova(gm1a)

## 5.4.2 Model predictions
data(mackp)
data(coast)
mackp$log.net.area <- 0*mackp$lon # make offset column
lon<-seq(-15,-1,1/4);lat<-seq(44,58,1/4)
zz<-array(NA,57*57)
zz[mackp$area.index]<-predict(gm1a,mackp)

image(lon,lat,matrix(zz,57,57),col=gray(0:32/32),
      cex.lab=1.5,cex.axis=1.4)
contour(lon,lat,matrix(zz,57,57),add=TRUE)
lines(coast$lon,coast$lat,col=1)

library(MASS)
br1 <- mvrnorm(n=1000,coef(gm1a),gm1a$Vp)
Xp <- predict(gm1a,newdata=mackp,type="lpmatrix")
mean.eggs1 <- colMeans(exp(Xp%*%t(br1)))
hist(mean.eggs1)

f<-fitted(gm1a)
form<-egg.count~offset(log.net.area)+s(lon,lat,bs="ts",k=100)+
                s(I(b.depth^.5),bs="ts")+s(temp.20m,bs="ts")
mack.bs <- mack
n <- nrow(mack)
br <- matrix(0,0,length(coef(gm1a)))
for (i in 1:19)
{ e <- rpois(rep(1,n),f) - f
  y <- round(f+e*gm1a$sig2^.5) ## deal with overdispersion
  y[y<0] <- 0
  mack.bs$egg.count <- y
  sp <-
  gam(form,data=mack.bs,family=poisson,scale=-1,gamma=1.4)$sp
  b <- gam(form,data=mack,family=poisson,sp=sp,scale=-1)
  br <- rbind(br,mvrnorm(n=100,coef(b),b$Vp))
}
br <- rbind(br,mvrnorm(n=100,coef(gm1a),gm1a$Vp))
mean.eggs <- colMeans(exp(Xp%*%t(br)))
hist(mean.eggs)

## 5.5 Portuguese larks example

data(bird)
ind <- sample(1:25100,2000,replace=FALSE)
m2 <- gam(crestlark ~ s(x,y,k=100),data=bird,family=binomial,
          knots=bird[ind,],gamma=1.4)
m2

bird$n <- bird$y*0+1 # when summed, gives binomial denominator
bird$n[is.na(bird$crestlark)] <- NA
ba <- aggregate(data.matrix(bird),by=list(bird$QUADRICULA),
                FUN=sum,na.rm=TRUE)
ba$x <- ba$x / 25   # don't want locations summed!
ba$y <- ba$y / 25

m10 <- gam(cbind(crestlark,n-crestlark)~s(x,y,k=100),
             data=ba,family=binomial,gamma=1.4)
m10

library(geoR)
coords<-matrix(0,1004,2);coords[,1]<-ba$x;coords[,2]<-ba$y
gb<-list(data=residuals(m10,type="d"),coords=coords)
plot(variog(gb,max.dist=1e5))
plot(fitted(m10),residuals(m10))

mx <- sort(unique(bird$x));my<-sort(unique(bird$y))
nx<-length(mx);ny<-length(my)
ixm <- 1:nx;names(ixm)<-mx
ix <- ixm[as.character(bird$x)]
iym <- 1:ny;names(iym)<-my
iy <- iym[as.character(bird$y)]
um <- matrix(NA,nx,ny)
fv10 <- predict(m10,bird,type="response")
my<-my/1000;mx <- mx/1000
um[ix+(iy-1)*nx] <- fv10
image(mx,my,matrix(um,nx,ny),xlab="km east",ylab="km north")
contour(mx,my,matrix(um,nx,ny),add=TRUE)

## 5.6.1 Package `gam' (note need to run chicago pre-proc first)

detach(package:mgcv)  
library(gam)
library(akima) # needed for plotting
bfm <- gam(death~s(time,df=140)+lo(o3,tmp,span=.1),
           family=poisson,control=gam.control(bf.maxit=150))
summary(bfm)
plot(bfm,se=T,rug=F,phi=30,theta=-30)
detach(package:gam)

## 5.6.2 Package `gss'

library(gss)
ssm <- gssanova1(death~time+o3*tmp,family="poisson",nbasis=200)

summary(ssm)

tp <- seq(min(time),max(time),length=500)
fvt <- predict(ssm,newdata=data.frame(time=tp,
               o3=rep(mean(o3),500),tmp=rep(mean(tmp),500)),
               include=list("time"))
plot(tp,fvt,type="l",xlab="time",ylab="time effect")

m <- 40
o3m <- seq(min(o3),max(o3),length=m)
tmpm <- seq(min(tmp),max(tmp),length=m)
tmpp <- rep(tmpm,rep(m,m))
o3p <- rep(o3m,m)
pd <- data.frame(time=rep(0,m*m),o3=o3p,tmp=tmpp)
fv <- predict(ssm,newdata=pd,include=list("o3","tmp","o3:tmp"))
library(mgcv)
ind <- exclude.too.far(o3p,tmpp,o3,tmp,dist=0.04)
fv[ind] <- 0
persp(o3m,tmpm,matrix(fv,m,m),phi=30,theta=-30,zlab="o3*tmp",
      xlab="o3",ylab="tmp")
