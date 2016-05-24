### R code from vignette source 'Introduction.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Introduction.Rnw:75-77
###################################################
options(width=80,useFancyQuotes="UTF-8")
require(rstpm2)


###################################################
### code chunk number 2: Introduction.Rnw:79-92
###################################################
data(popmort, package="rstpm2")
data(colon,package="rstpm2")
popmort2 <- transform(popmort,exitage=age,exityear=year,age=NULL,year=NULL)
colon2 <- within(colon, {
  status <- ifelse(surv_mm>120.5,1,status)
  tm <- pmin(surv_mm,120.5)/12
  exit <- dx+tm*365.25
  sex <- as.numeric(sex)
  exitage <- pmin(floor(age+tm),99)
  exityear <- floor(yydx+tm)
  ##year8594 <- (year8594=="Diagnosed 85-94")
})
colon2 <- merge(colon2,popmort2)


###################################################
### code chunk number 3: Introduction.Rnw:96-99
###################################################
fit0 <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
              data=colon2,
              bhazard=colon2$rate, df=5)


###################################################
### code chunk number 4: Introduction.Rnw:101-106
###################################################
summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate,
                     df=5,cure=TRUE))
predict(fit,head(colon2),se.fit=TRUE)


###################################################
### code chunk number 5: Introduction.Rnw:124-128
###################################################
newdata.eof <- data.frame(year8594 = unique(colon2$year8594),
                          tm=max(colon2$tm)) 
1-predict(fit0, newdata.eof, type="surv", se.fit=TRUE)
1-predict(fit, newdata.eof, type="surv", se.fit=TRUE)


###################################################
### code chunk number 6: Introduction.Rnw:132-143
###################################################
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 85-94"), ylim=0:1,
     xlab="Time since diagnosis (years)", ylab="Relative survival")
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 75-84"),
     add=TRUE,line.col="red",rug=FALSE)
plot(fit,newdata=data.frame(year8594 = "Diagnosed 85-94"),
     add=TRUE,ci=FALSE,lty=2,rug=FALSE)
plot(fit,newdata=data.frame(year8594="Diagnosed 75-84"),
     add=TRUE,rug=FALSE,line.col="red",ci=FALSE,lty=2)
legend("topright",c("85-94 without cure","75-84 without cure",
                    "85-94 with cure","75-84 with cure"),
       col=c(1,2,1,2), lty=c(1,1,2,2), bty="n")


###################################################
### code chunk number 7: Introduction.Rnw:150-165
###################################################
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 85-94"), 
     ylim=c(0,0.5), type="hazard",
     xlab="Time since diagnosis (years)",ylab="Excess hazard")
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 75-84"),
     type="hazard",
     add=TRUE,line.col="red",rug=FALSE)
plot(fit,newdata=data.frame(year8594 = "Diagnosed 85-94"),
     type="hazard",
     add=TRUE,ci=FALSE,lty=2,rug=FALSE)
plot(fit,newdata=data.frame(year8594="Diagnosed 75-84"),
     type="hazard",
     add=TRUE,rug=FALSE,line.col="red",ci=FALSE,lty=2)
legend("topright",c("85-94 without cure","75-84 without cure",
                    "85-94 with cure","75-84 with cure"),
       col=c(1,2,1,2), lty=c(1,1,2,2), bty="n")


