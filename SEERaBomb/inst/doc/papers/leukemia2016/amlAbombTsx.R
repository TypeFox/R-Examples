# amlAbombTsx.R (time since exposure plot)
# Uses data of Hsu et al. The Incidence of Leukemia ... among Atomic Bomb
# Survivors: 1950-2001. Radiat. Res. 179, 361-382 (2013).
(3.5/9)^0.5 # 0.623 Gy is the average whole body risk equivalent of a typical radiation therapy
rm(list=ls()) 
library(dplyr)
# library(SEERaBomb)
# mkAbomb() # if you haven't made abomb.db yet
db <- src_sqlite("~/data/abomb/abomb.db")
src_tbls(db)    #see what tables are in it
collect(tbl(db, sql("SELECT * from hemeDesc"))) 
d=collect(tbl(db, sql("SELECT * from heme"))) 
head(d)
sapply(d,class)
(d=d%>%select(city:un4gy,agexg:calg,py,agex:year,AML:AMLtot,D))
# d=round(d,2)
negs=d%>%filter(D<0) 
sum(negs$AMLtot) # 16 AMLs with unknown doses
sum(negs$py) # 200k py
d=d%>%filter(D>=0) 
sum(d$AMLtot) # 176 AMLs (= eligible AML cases in Table 1 of Hsu et al)
sum(d$py) #3,613,406 is pretty close to 3,613,404 in Table 2 of Hsu et al  

d$Dc=cut(d$D,breaks=c(0,0.005,0.1,0.2,0.5,1,2,10),include.lowest=TRUE,right=F)
head(d,2)
d%>%group_by(Dc)%>%summarize(PY=sum(py),AML=sum(AMLtot)) # reproduces PY in Hsu Table 2
d%>%group_by(calg)%>%summarize(PY=sum(py),AML=sum(AMLtot)) 
# merge calg=14 with 13
d$calg[d$calg==14]=13
d%>%group_by(calg)%>%summarize(PY=sum(py),AML=sum(AMLtot))#last group has more PY now
(years=d%>%group_by(calg)%>%summarize(aveY=weighted.mean(year,py))) 
(Year=round(years$aveY,1)) #time points in years
(Tsx=round(years$aveY,1)-1945.7) #time since exposure in years, for plotting
head(d)
d=d%>%mutate(aml=AMLtot,sv=D,tsx=age-agex,city=city-1,sex=sex-1)%>%select(-(AML:AMLtot),-year,-Dc,-D)
d=d%>%mutate(svc=cut(sv,c(-1,.01,.4,10),labels=c("low","med","high")),agec=cut(age,seq(0,110,10)))
head(d,3)
agem=seq(5,110,10) # midpoints of age intervals

# pool across city and sex to get background expected values
(BKPY=with(d,tapply(py,list(sv=svc,age=agec),sum,na.rm=TRUE)))  
sum(BKPY)
bkCounts=with(d,tapply(aml,list(sv=svc,age=agec),sum,na.rm=TRUE)) 
(X=data.frame(age=agem,cnts=bkCounts[1,],py=BKPY[1,]))
library(mgcv)
bkg=gam(cnts ~ s(age)+offset(log(py)),family=poisson(),data=X,method="REML") 
X$Ecases=exp(predict(bkg))
X$Eincid=X$Ecases/X$py
X$incid=X$cnts/X$py
with(X,plot(age,Eincid,type="l"))
with(X,points(age,incid))

d$bk=as.numeric(exp(predict(bkg,newdata=d))/d$py) # numerator is expected cases
with(d,plot(age,bk))
X0=c(rep(-10,13))
flin<-function(L)  {
  d$mn = d$bk*(1+d$sv^2*exp(L[d$calg]))*d$py # relative risk quad in dose
  -sum(d$aml*log(d$mn) - d$mn,na.rm=T)
}
range(d$calg)
head(d)
sol=optim(X0,flin,method="L-BFGS-B",control=list(maxit=400));sol  # the idea here is to let the data speak through L, a bit like a one way anova
wait=exp(sol$par)+1
plot(Year,wait,pch=5,col="red",type="b",ylab="w(t)",main="IR-to-AML waiting time pdf shape")

# repeat with bbmle
library(bbmle)
d=transform(d,calg=as.factor(calg)) 
(s=summary(m<-mle2(aml~dpois(lambda=py*bk*(1+sv^2*exp(f))),
                   parameters=list(f~-1 + calg),
                   start=list(f=-5),data=d) )  )
ci=confint(m)
wait=exp(coef(s)[,1])
waitL=exp(coef(s)[,1]-1.96*coef(s)[,2])
waitU=exp(coef(s)[,1]+1.96*coef(s)[,2])
waitLp=exp(ci[,1])
waitUp=exp(ci[,2]) # profile versions
(dfc=data.frame(Tsx,data.frame(wait,waitL,waitU,waitLp,waitUp)+1 ))
dfc[is.na(dfc)]=0 # minus infinity maps to zero, so set them there

library(ggplot2)
quartz(width=6.5,height=4)
theme_update(plot.title = element_text(size = rel(2.3)),
             axis.title = element_text(size = rel(2.3)),
             axis.text = element_text(size = rel(2.3))) 
p=ggplot(dfc, aes(x=Tsx, y=wait)) +  geom_point(size=6) 
p=p+labs(x="Years since exposure",y="AML RR (1 Sv)") + 
  scale_x_continuous(breaks=seq(0,50,by=10),limits=c(0,55))
p=p+geom_errorbar(aes(ymin=waitLp, ymax=waitUp),width=.01)
p+geom_abline(intercept=1, slope=0)
ggsave("~/Results/amlMDS/RRtsxAMLabomb.png")

