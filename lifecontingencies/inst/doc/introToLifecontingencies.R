## ----setup, include=FALSE------------------------------------------------
library(knitr)
rm(list=ls())

## ----load, echo=TRUE, warning=FALSE, message=FALSE-----------------------
library(lifecontingencies) #load the package

## ----finmat1, echo=TRUE, warning=FALSE, message=FALSE--------------------
capitals <- c(-1000,200,500,700)
times <- c(0,1,2,5)
#calculate a present value
presentValue(cashFlows=capitals, timeIds=times, 
             interestRates=0.03)

## ----finmat2, echo=TRUE, warning=FALSE, message=FALSE--------------------
ann1 <- annuity(i=0.03, n=5, k=1, type="immediate")
ann2 <- annuity(i=0.03, n=5, k=12, type="due")
c(ann1,ann2)

## ----finmat3, echo=TRUE, warning=FALSE, message=FALSE--------------------
bondPrice<-5*annuity(i=0.03,n=10)+100*1.03^-10
bondPrice

## ----demo1---------------------------------------------------------------
#create an demo lifetable
xDemo<-seq(from=0,to=5,by=1)
lxDemo<-c(100,95,90,60,30,5)
lifetableDemo<-new("lifetable",x=xDemo,
                   lx=lxDemo,name="Demo")

## ----demo2---------------------------------------------------------------
data(demoIta) #using the internal Italian LT data set
lxIPS55M <- with(demoIta, IPS55M)
#performing some fixings
pos2Remove <- which(lxIPS55M %in% c(0,NA))
lxIPS55M <-lxIPS55M[-pos2Remove]
xIPS55M <-seq(0,length(lxIPS55M)-1,1)
#creating the table
ips55M <- new("lifetable",x=xIPS55M, 
lx=lxIPS55M,name="IPS 55 Males")

## ----demo3, tidy=TRUE----------------------------------------------------
#decrements between age 65 and 70
dxt(ips55M, x=65, t = 5)
#probabilities of death between age 80 and 85
qxt(ips55M, x=80, t=2)
#expected curtate lifetime
exn(ips55M, x=65) 

## ----createacttable, tidy=FALSE------------------------------------------
#creates a new actuarial table
ips55Act<-new("actuarialtable",
x=ips55M@x,lx=ips55M@lx,
interest=0.02,name="IPS55M")

## ----pureend, tidy=FALSE-------------------------------------------------
#compute APV
APV=50e3*Exn(actuarialtable = 
ips55Act,x=30,n=35)
#compute Premium
P=APV/axn(actuarialtable =
ips55Act,x=30,n=35)
c(APV,P)

## ----endowmentcalcs, tidy=FALSE------------------------------------------
#defining the ranges
interest.range<-seq(from=0.015, to=0.035,by=0.001)
term.range<-seq(from=20, to=40,by=1)
#computing APV sensitivities
apv.interest.sensitivity<-sapply(interest.range,
FUN = "AExn",actuarialtable=ips55Act,x=30,n=30)
apv.term.sensitivity<-sapply(term.range,FUN = "AExn",
                             actuarialtable=ips55Act,x=30)

## ----endowmentplot, tidy=FALSE, echo=FALSE,fig.width=5, fig.height=5,fig.align='center'----
par(mfrow=c(1,2))
plot(x=interest.range, y=apv.interest.sensitivity,type="l",xlab="interest rate",ylab="APV",main="APV by Interest Rate")
plot(x=term.range, y=apv.term.sensitivity,type="l",xlab="term",ylab="APV",main="APV by term")

## ----reserves, tidy=FALSE------------------------------------------------
#compute the APV and premium
APV=100e3*Axn(actuarialtable = ips55Act,x=25,n=40) 
P=APV/axn(actuarialtable = ips55Act,x=25,n=40)
#define a reserve function
reserveFunction<-function(t) 
  100e3*Axn(actuarialtable = ips55Act,x=25+t,n=40-t) - 
  P *axn(actuarialtable = ips55Act,x=25+t,n=40-t)
reserve<-sapply(0:40,reserveFunction)

## ----reserves2, tidy=FALSE, echo=FALSE, fig.align='center'---------------
plot(x=0:40,y=reserve,main="Reserve",
     xlab="Policy Age",ylab="Reserve outstanding",type="l")

## ----AEXn1 , tidy=FALSE, size="small"------------------------------------
#analyzing and Endowment of 100K on x=40, n=25
#compute APV
APV=AExn(actuarialtable = ips55Act,x=40,n=25) 
#sampling
AEXnDistr<-rLifeContingencies(n=10e3,
lifecontingency = "AExn",x = 40,
t=25,object = ips55Act)

## ----AExn1, tidy=FALSE, size="small"-------------------------------------
#assess if the expected value match the theoretical one
t.test(x=AEXnDistr,mu = APV)

## ----AEXn2, tidy=FALSE, echo=FALSE---------------------------------------
hist(AEXnDistr, main="Endowment Actuarial Value Distribution",
     probability = TRUE, col="steelblue")

## ----leecarter01, tidy=FALSE, include=FALSE, results='hide'--------------
#library(demography)
#italy.demo<-hmd.mx("ITA", username="spedicato_giorgio@yahoo.it", password="mortality")

## ----leecarter0, tidy=FALSE, warning=FALSE, message=FALSE----------------
#load the package and the italian tables
library(demography) 
#italyDemo<-hmd.mx("ITA", username="yourUN", 
#password="yourPW")
load(file="mortalityDatasets.RData") #load the dataset

## ----leecarter1, tidy=FALSE, warning=FALSE-------------------------------
#calibrate lee carter
italy.leecarter<-lca(data=italyDemo,series="total",
                     max.age=103,adjust = "none")
#perform modeling of kt series
kt.model<-auto.arima(italy.leecarter$kt)
#projecting the kt
kt.forecast<-forecast(kt.model,h=100) 

## ----leecarter2, tidy=FALSE, size='tiny'---------------------------------
#indexing the kt
kt.full<-ts(union(italy.leecarter$kt, kt.forecast$mean),
            start=1872)  
#getting and defining the life tables matrix
mortalityTable<-exp(italy.leecarter$ax
+italy.leecarter$bx%*%t(kt.full)) 
rownames(mortalityTable)<-seq(from=0, to=103)
colnames(mortalityTable)<-seq(from=1872, 
to=1872+dim(mortalityTable)[2]-1)

## ----leecarter2plot, tidy=FALSE, echo=FALSE------------------------------
plot.ts(kt.full, main="historical and projected KT",xlab="year",
        ylab="kt",col="steelblue")
abline(v=2009,col="darkred",lwd=2.5)

## ----leecarter3, tidy=FALSE----------------------------------------------
getCohortQx<-function(yearOfBirth)
{
  colIndex<-which(colnames(mortalityTable)
                  ==yearOfBirth) #identify 
  #the column corresponding to the cohort 
  #definex the probabilities from which 
  #the projection is to be taken
  maxLength<-min(nrow(mortalityTable)-1,
                 ncol(mortalityTable)-colIndex)
  qxOut<-numeric(maxLength+1)
  for(i in 0:maxLength)
    qxOut[i+1]<-mortalityTable[i+1,colIndex+i]
  #fix: we add a fictional omega age where 
  #death probability = 1
  qxOut<-c(qxOut,1)
  return(qxOut)
}

## ----leecarter4, tidy=FALSE, size='scriptsize'---------------------------
#generate the life tables
qx1920<-getCohortQx(yearOfBirth = 1920)
lt1920<-probs2lifetable(probs=qx1920,type="qx",
name="Table 1920")
at1920<-new("actuarialtable",x=lt1920@x,
lx=lt1920@lx,interest=0.015)
qx1950<-getCohortQx(yearOfBirth = 1950)
lt1950<-probs2lifetable(probs=qx1950,
type="qx",name="Table 1950")
at1950<-new("actuarialtable",x=lt1950@x,
lx=lt1950@lx,interest=0.015)
qx1980<-getCohortQx(yearOfBirth = 1980)
lt1980<-probs2lifetable(probs=qx1980,
type="qx",name="Table 1980")
at1980<-new("actuarialtable",x=lt1980@x,
lx=lt1980@lx,interest=0.015)

## ----leecarter5, tidy=FALSE, echo=TRUE-----------------------------------
cat("Results for 1920 cohort","\n")
c(exn(at1920,x=65),axn(at1920,x=65))
cat("Results for 1950 cohort","\n")
c(exn(at1950,x=65),axn(at1950,x=65))
cat("Results for 1980 cohort","\n")
c(exn(at1980,x=65),axn(at1980,x=65))

