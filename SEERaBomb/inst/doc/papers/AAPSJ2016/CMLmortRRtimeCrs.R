library(ggplot2)
library(dplyr)
load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
d=canc%>%filter(cancer=="CML")
d=d%>%mutate(status=as.numeric(COD>0),surv=(surv+0.5)/12)%>%filter(surv<200) #remove survival unknowns =9999 months 
# d=d%>%mutate(Period=cut(yrdx,breaks=c(1972,1985,1995,2000,2005,2012),dig.lab=4))
d=d%>%mutate(Period=cut(yrdx,breaks=c(1972,1990,2002,2013),dig.lab=4))
d=d%>%select(Period,yrdx,agedx,sex,surv,status)
head(d)
L=split(d,d$Period)
str(L)
# library(demography)
# d=hmd.mx("USA", "username", "password") #make an account and put your info in here
# mrt=d$rate
# save(mrt,file="~/data/usMort/mrt.RData")
load("~/data/usMort/mrt.RData")
d%>%group_by(sex)%>%summarize(med=median(agedx))
mean(mrt$female["62",paste(1973:2013)])
mean(mrt$male["60",paste(1973:2013)])
2.5*1.5
2.5*1

brks=c(0,0.5,1,2,3,4,5,6,8) 
library(SEERaBomb)
LO=lapply(L,function (x) msd(x,mrt,brks))
str(LO)
for (i in names(LO)) LO[[i]]$Period=i
D=NULL
for (i in names(LO)) D=rbind(D,LO[[i]])
D

quartz(width=6,height=5)
theme_update(legend.position = c(.88, .39),  
             # legend.position="top",
             axis.text=element_text(size=rel(1.4)),
             axis.title=element_text(size=rel(1.4)),
             legend.title=element_text(size=rel(.9)),
             legend.text=element_text(size=rel(.9)),
             strip.text = element_text(size = rel(1.5)))
g=qplot(x=t,y=RR,data=D,col=Period,geom=c("line","point"), ylim=c(0,18),
        xlab="Years Since CML Diagnosis",ylab="Relative Risk of Mortality")
g=g+geom_line(size=1)+facet_grid(sex~.)+geom_abline(intercept=1, slope=0) 
g+geom_errorbar(aes(ymin=rrL,ymax=rrU,width=.15)) 
ggsave("~/Results/CML/mortRRtimeCrsTrends.png")  

  
             