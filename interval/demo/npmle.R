# compare to Icens package which has several different methods for calculating the npmle
library(interval)
data(bcos)
library(Icens)

# check either default or formula methods
#fitInt<-icfit(bcos$left,bcos$right)
ic1<-icfit(Surv(left,right,type="interval2")~1,data=bcos,control=icfitControl(epsilon=10^-12))
EM1<-EM(bcos[,c("left","right")],maxiter=10^4,tol=10^-12)
## all npmle estimates match to within 4 decimal places
check.fits<-function(fitic,fitEM,digits=5){
    class(fitEM)<-"icfit"
    em<-summary(fitEM,digits=digits)
    sEM<-em[em[,"Probability"]>10^-digits,]
    sic<-summary(fitic,digits=digits)
    sic<-sic[sic[,"Probability"]>10^-digits,]
    if (dim(sEM)[[1]]==dim(sic)[[1]]) out<- all(sEM[,"Probability"]==sic[,"Probability"])
    else out<-FALSE
    names(out)<-paste("Probs Match to ",digits,"decimal places")
    out 
}
# match to 8 decimal places
check.fits(ic1,EM1,digits=8)

# Simulate different data
set.seed(1232)
n<-10
nsim<-30
m<-rep(NA,nsim)
for (i in 1:nsim){
    left<-rpois(n,5)
    right<-left+1+rpois(n,2)
    if (i==21) dat<-data.frame(left,right)
    fitEM<-EM(cbind(left,right),maxit=10^4,tol=10^-12)
    fitic<-icfit(left,right,control=icfitControl(epsilon=10^-12,maxit=10^4))
    m[i]<-check.fits(fitic,fitEM,digits=8)
}
# note that the 21st data set does not match, saved it
(1:30)[!m]



# need to increase maxit for the 21st data set
# took too long so commented it out
#fitEM<-EM(dat[,c("left","right")],maxit=10^6,tol=10^-12)
#fitic<-icfit(dat$left,dat$right,control=icfitControl(epsilon=10^-12,maxit=10^6))

# still only match to 5 decimal places...sometimes EM is very slow
#check.fits(fitic,fitEM,digits=5)
#> check.fits(fitic,fitEM,digits=5)
#  (Left Right] Probability
#1     4      5        0.50
#2     5      6        0.00
#3     6      7        0.25
#4     7      9        0.25
#  (Left Right] Probability
#1     4      5        0.50
#2     5      6        0.00
#3     6      7        0.25
#4     7      9        0.25
#Probs Match to  5 decimal places 
#                            TRUE 
#Warning messages:
#1: In summary.icfit(fitEM, digits = digits) : fit did not converge
#2: In summary.icfit(fitic, digits = digits) : fit did not converge

#check.fits(fitic,fitEM,digits=6)
#> check.fits(fitic,fitEM,digits=6)
#  (Left Right] Probability
#1     4      5    0.499999
#2     5      6    0.000002
#3     6      7    0.249999
#4     7      9    0.250000
#  (Left Right] Probability
#1     4      5    0.499998
#2     5      6    0.000002
#3     6      7    0.249999
#4     7      9    0.250000
#Probs Match to  6 decimal places 
#                           FALSE 
#Warning messages:
#1: In summary.icfit(fitEM, digits = digits) : fit did not converge
#2: In summary.icfit(fitic, digits = digits) : fit did not converge
#> 

