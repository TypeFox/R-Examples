### R code from vignette source 'semiparametric-unemployment.Rnw'

###################################################
### code chunk number 1: semiparametric-unemployment.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: semiparametric-unemployment.Rnw:19-22 (eval = FALSE)
###################################################
## library(catdata)
## data(unemployment, package="catdata")
## attach(unemployment)


###################################################
### code chunk number 3: semiparametric-unemployment.Rnw:27-28 (eval = FALSE)
###################################################
## library(mgcv)


###################################################
### code chunk number 4: semiparametric-unemployment.Rnw:33-36 (eval = FALSE)
###################################################
## durbin[durbin==2] <- 0
## 
## gamage <- gam(durbin ~ s(age, bs="ps", m=c(2,1), k=15), family=binomial())


###################################################
### code chunk number 5: semiparametric-unemployment.Rnw:41-46 (eval = FALSE)
###################################################
## minage <- min(age)
## maxage <- max(age)
## ageindex <- seq(from=minage, to=maxage, by=0.01)
## 
## pred <- predict(gamage, newdata=data.frame("age"=ageindex), type="response")


###################################################
### code chunk number 6: semiparametric-unemployment.Rnw:51-59 (eval = FALSE)
###################################################
## bspline<-function(x,k,i,m=2)
## {if (m==-1)
## {res<-as.numeric(x<k[i+1]&x>=k[i])}
## else{
## z0<-(x-k[i])/(k[i+m+1]-k[i])
## z1<-(k[i+m+2]-x)/(k[i+m+2]-k[i+1])
## res<- z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1)}
## res}


###################################################
### code chunk number 7: semiparametric-unemployment.Rnw:64-74 (eval = FALSE)
###################################################
## k <- gamage$smooth[[1]]$knots
## 
## meanage <- c()
## 
## for (i in minage:maxage){
## meanage[i] <- sum(durbin[age==i])
## if(sum(durbin[age==i])!=0){
## meanage[i] <- mean(durbin[age==i])
## }
## }


###################################################
### code chunk number 8: semiparametric-unemployment.Rnw:79-86 (eval = FALSE)
###################################################
## par(cex=1.3, lwd=1.5)
## plot(ageindex, pred, type="l", ylim=c(0,0.8), col="gray", xlab="age/years", 
##      ylab="short-term unemployment")
## for(i in 1:15)
## {lines(ageindex,(0.5*gamage$coefficients[i+1]+0.7)*(bspline(x=ageindex,k=k,
## i=i+1,m=2)),type="l", col="gray")}
## points(minage:maxage, meanage[minage:maxage], pch=20)


###################################################
### code chunk number 9: semiparametric-unemployment.Rnw:91-94 (eval = FALSE)
###################################################
## gamage2 <- gam(durbin ~ s(age, bs="ps", fx=TRUE, m=c(2,1),k=15), family=binomial())
## 
## pred2 <- predict(gamage2, newdata=data.frame("age"=ageindex), type="response")


###################################################
### code chunk number 10: semiparametric-unemployment.Rnw:99-106 (eval = FALSE)
###################################################
## par(cex=1.3, lwd=1.5)
## plot(ageindex, pred2, type="l", ylim=c(0,0.8), col="gray", xlab="age/years", 
##      ylab="short-term unemployment")
## for(i in 1:15)
## {lines(ageindex, (0.1*gamage2$coefficients[i+1]+0.3)*
##   bspline(x=ageindex,k=k,i=i+1,m=2),type="l", col="gray")}
## points(minage:maxage, meanage[minage:maxage],pch=20)


###################################################
### code chunk number 11: semiparametric-unemployment.Rnw:109-111 (eval = FALSE)
###################################################
## detach(unemployment)
## rm(unemployment)


