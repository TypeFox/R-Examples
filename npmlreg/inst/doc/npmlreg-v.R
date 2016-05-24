### R code from vignette source 'npmlreg-v.Rtex'

###################################################
### code chunk number 1: npmlreg-v.Rtex:176-180
###################################################
data(galaxies, package="MASS")
galaxies[78]<-26960
gal<-as.data.frame(galaxies)
rm(galaxies)


###################################################
### code chunk number 2: npmlreg-v.Rtex:187-189
###################################################
gal$v1000<- gal$galaxies/1000
gal$v1000


###################################################
### code chunk number 3: npmlreg-v.Rtex:195-196
###################################################
library(npmlreg)


###################################################
### code chunk number 4: npmlreg-v.Rtex:201-202
###################################################
glm(v1000~1,data=gal)


###################################################
### code chunk number 5: npmlreg-v.Rtex:206-207
###################################################
(galaxy.np1 <- alldist(v1000~1,random=~1,random.distribution='np',k=1,data=gal))


###################################################
### code chunk number 6: npmlreg-v.Rtex:221-222
###################################################
galaxy.np1$dev


###################################################
### code chunk number 7: npmlreg-v.Rtex:229-232
###################################################
(galaxy.np2 <- alldist(v1000~1,random=~1,random.distribution='np',k=2,data=gal))
(galaxy.np3 <- alldist(v1000~1,random=~1,random.distribution='np',k=3,data=gal))
(galaxy.np4 <- alldist(v1000~1,random=~1,random.distribution='np',k=4,data=gal))


###################################################
### code chunk number 8: npmlreg-v.Rtex:242-243
###################################################
plot(galaxy.np4, plot.opt=3)


###################################################
### code chunk number 9: npmlreg-v.Rtex:251-256
###################################################
(galaxy.np5 <- alldist(v1000~1,random=~1,k=5,data=gal, verbose=FALSE))$disp
(galaxy.np6 <- alldist(v1000~1,random=~1,k=6,tol=0.2,data=gal,verbose=FALSE))$disp
(galaxy.np7 <- alldist(v1000~1,random=~1,k=7,tol=0.12,data=gal,verbose=FALSE))$disp
(galaxy.np8 <- alldist(v1000~1,random=~1,k=8,tol=0.2,data=gal,verbose=FALSE))$disp
(galaxy.np9 <- alldist(v1000~1,random=~1,k=9,tol=0.06,data=gal,verbose=FALSE))$disp


###################################################
### code chunk number 10: npmlreg-v.Rtex:282-283
###################################################
summary(galaxy.np4u <- alldist(v1000~1, random=~1, k=4, tol=0.5, data=gal, lambda=1, verbose=FALSE))


###################################################
### code chunk number 11: npmlreg-v.Rtex:289-290 (eval = FALSE)
###################################################
## plot(galaxy.np4u, plot.opt=15, height=5)


###################################################
### code chunk number 12: npmlreg-v.Rtex:306-307
###################################################
plot(galaxy.np4u, plot.opt=15)


###################################################
### code chunk number 13: npmlreg-v.Rtex:338-340
###################################################
(galaxy.np8us <- alldist(v1000~1, random=~1, k=8, tol=0.5, data=gal, lambda=1, verbose=FALSE, spike.protect=TRUE))
galaxy.np8us$sdev$sdevk


###################################################
### code chunk number 14: npmlreg-v.Rtex:351-353
###################################################
(galaxy.np8ud <- alldist(v1000~1, random=~1, k=8, tol=0.5, data=gal, lambda=0.99))
galaxy.np8ud$sdev$sdevk


###################################################
### code chunk number 15: npmlreg-v.Rtex:364-366
###################################################
par(mfrow=c(1,1), cex=0.65)
tolfind(v1000~1, random=~1, k=8, data=gal, lambda=1, find.in.range=c(0.0,0.6), steps=12, plot.opt=0, verbose=FALSE, noformat=TRUE)[c(3,4)]


###################################################
### code chunk number 16: npmlreg-v.Rtex:379-382
###################################################
data(fabric)
(faults0 <- glm(y ~ 1, family=poisson(link=log),data=fabric))
(faults1 <- glm(y ~ x, family=poisson(link=log),data=fabric)) 


###################################################
### code chunk number 17: npmlreg-v.Rtex:401-404
###################################################
(faults.g1<- alldist(y ~ x, family=poisson(link=log), random=~1, data= fabric,k=1, random.distribution="gq")) 
(faults.g2<- alldist(y ~ x, family=poisson(link=log), random=~1, data= fabric,k=2, random.distribution="gq")) 
(faults.g3<- alldist(y ~ x, family=poisson(link=log), random=~1, data= fabric,k=3, random.distribution="gq",verbose=F)) 


###################################################
### code chunk number 18: npmlreg-v.Rtex:408-409
###################################################
faults.g1$dev


###################################################
### code chunk number 19: npmlreg-v.Rtex:417-419
###################################################
(faults.np2<- alldist(y ~ x, family=poisson(link=log), random=~1, data= fabric,k=2, random.distribution="np")) 
(faults.np3<- alldist(y ~ x, family=poisson(link=log), random=~1, data= fabric,k=3, random.distribution="np",verbose=FALSE)) 


###################################################
### code chunk number 20: npmlreg-v.Rtex:428-429
###################################################
predict(faults.g2, type="response",newdata=fabric[1:6,])


###################################################
### code chunk number 21: npmlreg-v.Rtex:433-434
###################################################
predict(faults.g2, type="response")[1:6]


###################################################
### code chunk number 22: npmlreg-v.Rtex:447-450
###################################################
data(rainfall, package="forward")
rainfall$x<-rainfall$Rain/1000  
rainfall$x2<- rainfall$x^2; rainfall$x3<- rainfall$x^3


###################################################
### code chunk number 23: npmlreg-v.Rtex:457-458
###################################################
 (toxo.np<- alldist(cbind(Cases,Total-Cases)~1, random=~1, data=rainfall, k=3, family=binomial(link=logit)))


###################################################
### code chunk number 24: npmlreg-v.Rtex:462-463
###################################################
 toxo.np$disparity


###################################################
### code chunk number 25: npmlreg-v.Rtex:470-471
###################################################
 (toxo.npx<- alldist(cbind(Cases,Total-Cases)~x, random=~1, data=rainfall, k=3, family=binomial(link=logit)))


###################################################
### code chunk number 26: npmlreg-v.Rtex:476-477
###################################################
 (toxo.npxx<- alldist(cbind(Cases,Total-Cases)~x, random=~x, data=rainfall, k=3, family=binomial(link=logit)))


###################################################
### code chunk number 27: npmlreg-v.Rtex:484-485
###################################################
  round(t(toxo.np$post.prob),digits=2)


###################################################
### code chunk number 28: npmlreg-v.Rtex:495-497
###################################################
round(toxo.ebp<-toxo.np$ebp,digits=3)
round(exp(toxo.ebp)/(1+exp(toxo.ebp)),digits=4)


###################################################
### code chunk number 29: npmlreg-v.Rtex:501-502
###################################################
predict(toxo.np, type="response")


###################################################
### code chunk number 30: npmlreg-v.Rtex:505-506
###################################################
fitted(toxo.np)


###################################################
### code chunk number 31: npmlreg-v.Rtex:511-512
###################################################
predict(toxo.npx,type="response",newdata=data.frame(x=2))


###################################################
### code chunk number 32: npmlreg-v.Rtex:526-528
###################################################
data(hosp)
(fitnp3<-  alldist(duration~age+temp1, data=hosp,k=3, family=Gamma(link=log),tol=0.2)) 


###################################################
### code chunk number 33: npmlreg-v.Rtex:532-533
###################################################
 fitnp3$shape


###################################################
### code chunk number 34: npmlreg-v.Rtex:539-540
###################################################
 (fitnp3e<-  alldist(duration~age+temp1, data=hosp,k=3, family=Gamma(link=log),tol=0.2,shape=1))


###################################################
### code chunk number 35: npmlreg-v.Rtex:579-583
###################################################
 data(Oxboys, package = "nlme")
 Oxboys$boy <- gl(26,9) 
 plot(Oxboys$age[Oxboys$boy==1],Oxboys$height[Oxboys$boy==1],ylim=c(125,175),type='b',pch=1,xlab='age',ylab='height')
 for (i in 2:nlevels(Oxboys$Subject)){lines(Oxboys$age[Oxboys$boy==i],Oxboys$height[Oxboys$boy==i], pch=1,type='b',col=i)}


###################################################
### code chunk number 36: npmlreg-v.Rtex:595-596
###################################################
 (Oxboys.g20 <- allvc(height~age,random=~1|boy,data=Oxboys,random.distribution='gq',k=20))


###################################################
### code chunk number 37: npmlreg-v.Rtex:604-605
###################################################
Oxboys.g20$rsdev^2/(Oxboys.g20$rsdev^2+ Oxboys.g20$sdev$sdev^2)


###################################################
### code chunk number 38: npmlreg-v.Rtex:609-611
###################################################
 (Oxboys.np7 <- allvc(height~age,random=~1|boy,data=Oxboys,random.distribution='np',k=7))
 (Oxboys.np8 <- allvc(height~age,random=~1|boy,data=Oxboys,random.distribution='np',k=8)) 


###################################################
### code chunk number 39: npmlreg-v.Rtex:619-620
###################################################
 plot(Oxboys.np8, plot.opt=2)


###################################################
### code chunk number 40: npmlreg-v.Rtex:626-627
###################################################
 (Oxboys.np8s <- allvc(height~age,random=~age|boy,data=Oxboys,random.distribution='np',k=8))


###################################################
### code chunk number 41: npmlreg-v.Rtex:631-632
###################################################
  Oxboys.np8$disp-Oxboys.np8s$disp


###################################################
### code chunk number 42: npmlreg-v.Rtex:639-640
###################################################
 plot(Oxboys.np8, plot.opt=2)


###################################################
### code chunk number 43: npmlreg-v.Rtex:654-655
###################################################
 data(irlsuicide)


###################################################
### code chunk number 44: npmlreg-v.Rtex:673-674
###################################################
citation(package="npmlreg")


###################################################
### code chunk number 45: npmlreg-v.Rtex:753-754
###################################################
ls("package:npmlreg")


