runall<-FALSE
if(runall==TRUE){
#load("../Data/BlackwellData")

blackwell<-read.table("../Data/Blackwell.tab",header=TRUE)



##New help file
form1<-"d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 + camp.length + camp.length + deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office"

fit1<-CBMSM(formula = form1, time=blackwell$time,id=blackwell$demName,data=blackwell, type="MSM",  iterations = NULL, twostep = TRUE, msm.variance = "approx", time.vary = TRUE)

bal1<-balance.CBMSM(fit1)

##Effect estimation

lm(demprcnt[time==1]~fit1$treat.hist,data=blackwell,weights=fit1$weights)
lm(demprcnt[time==1]~fit1$treat.cum,data=blackwell)





##Attempting to replicate



#write.table(file="../Data/Blackwell.tab",blackwell)
blackwell<-read.table("../Data/Blackwell.tab")
form1<-"d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 + camp.length + camp.length + deminc + base.poll + as.factor(year) +   base.und + office"

fit1<-CBMSM(formula = form1, time=time,id=id,data=blackwell, type="MSM",  iterations = NULL, twostep = TRUE, msm.variance = "full", time.vary = TRUE)

fit2<-CBMSM(formula = form1, time=time,id=id,data=blackwell, type="MSM",  iterations = NULL, twostep = TRUE, msm.variance = "approx", time.vary = TRUE)

dv<-blackwell$demprcnt[blackwell$time==1]
treat.mat<-sapply(1:5,FUN=function(x) treat[time==x])
treat.cum<-rowSums(treat.mat)
colnames(treat.mat)<-paste("treat_",1:5,sep="")

lm(dv~treat.mat,w=fit1$weights[blackwell$time==1])
lm(dv~treat.mat,w=fit2$weights[blackwell$time==1])

lm(dv~treat.cum,w=fit1$weights[time==1])
lm(dv~treat.cum,w=fit2$weights[time==1])



}
