# prosBreastTimeCrs.R   (RR time courses after Prostate and Breast first cancers)
rm(list=ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(SEERaBomb)
system.time(load("~/Results/CLL/pm.RData")) # 4 secs to load. 
system.time(load("~/Results/CLL/pf.RData")) # 4 secs to load. 
Dm=mkDF(pm)
Df=mkDF(pf)
Dm=Dm%>%filter(cancer1=="prostate")
Df=Df%>%filter(cancer1=="breast")
Dm$cancer="Prostate"
Df$cancer="Breast"
D=rbind(Dm,Df)
D$Radiation="No"
D$Radiation[D$trt=="rad"]="Yes"
D$Period="Early"
D$Period[D$t>5]="Late"
D=D%>%select(int,t,RR,L=rrL,U=rrU,Radiation,Period,cancer)
head(D)

graphics.off()
quartz(width=7,height=3.7)

theme_set(theme_bw())
theme_update(legend.position = c(.65, .45),
             axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.title=element_text(size=rel(0.9)),
             legend.text=element_text(size=rel(0.9)),
             strip.text = element_text(size = rel(1.5)))
g=qplot(x=t,y=RR,data=D,col=Radiation,geom=c("line","point"),#xlim=c(-.1,24),
        xlab="Years Since Diagnosis of First Cancer",ylab="CLL Relative Risk")
g=g+scale_y_log10(breaks=c(0.3,1,10),labels=c("0.3","1","10"),limits=c(0.15,30))
g=g+facet_grid(cancer~Period,scales="free")+geom_abline(intercept=0, slope=0)
g=g+  geom_errorbar(aes(ymin=L,ymax=U,width=.05))
g = g + scale_color_grey(start = 0, end = 0.6)
g
ggsave("~/Results/CLL/prosBrsTimCrs.eps")  
