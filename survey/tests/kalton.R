library(survey)

ab<-expand.grid(a=factor(1:4),b=factor(1:3))

kaltonsample<-ab[rep(1:12,c(20,50,100,30,40,140,50,100,40,310,50,70)),]

kaltonpop<-ab[rep(1:12,c(80,60,170,55,40,150,60,165,55,340,200,125)),]

jointpop<-colSums(model.matrix(~a*b,kaltonpop))
marginalpop<-colSums(model.matrix(~a+b,kaltonpop))
gregpop<-colSums(model.matrix(~as.numeric(a)+as.numeric(b),kaltonpop))

dkalton<-svydesign(id=~1,data=kaltonsample)

dps<-postStratify(dkalton,~a+b,xtabs(~a+b,kaltonpop))

drake<-rake(dkalton, list(~a,~b),list(xtabs(~a,kaltonpop),xtabs(~b,kaltonpop)),control=list(epsilon=0.0001))

dcalps<-calibrate(dkalton, ~a*b, jointpop)
dcalrake<-calibrate(dkalton,~a+b, marginalpop, calfun="raking")
dlinear<-calibrate(dkalton, ~a+b, marginalpop)

dtrunclinear<-calibrate(dkalton, ~a+b, marginalpop,bounds=c(0.5,2.2))

dlogit<-calibrate(dkalton, ~a+b, marginalpop,bounds=c(0.5,2.2),calfun="logit")

dgreg<-calibrate(dkalton,~as.numeric(a)+as.numeric(b), gregpop)


#table A
 round(svytable(~a+b,dps)/xtabs(~a+b,kaltonsample),2)
 round(svytable(~a+b,dcalps)/xtabs(~a+b,kaltonsample),2)

#table B
 round(svytable(~a+b,drake)/xtabs(~a+b,kaltonsample),2)
 round(svytable(~a+b,dcalrake)/xtabs(~a+b,kaltonsample),2)

#table C
round(svytable(~a+b,dlinear)/xtabs(~a+b,kaltonsample),2)

#table D
round(svytable(~a+b,dgreg)/xtabs(~a+b,kaltonsample),2)

#table G
round(svytable(~a+b,dlogit)/xtabs(~a+b,kaltonsample),2)

#table G
round(svytable(~a+b,dtrunclinear)/xtabs(~a+b,kaltonsample),2)
