# timeCrs.R   (CLL RR time courses after Non-Hematological first cancers)
rm(list=ls()) 
library(SEERaBomb) 
library(reshape2)
library(dplyr)  
library(grid)
load("~/Results/CLL/pm.RData") 
load("~/Results/CLL/pf.RData") 
brks=c(0,.1,.2,.3,0.6,1,1.5,2,2.5,3,4,5,7,10,13,16,20) # res back up since now only CLL
(brkS=paste0("b",paste(brks,collapse="_")))
HM=c("AML","AMLti","APL","MDS","CMML","CML","MPN","ALL","CLL","HCL","OL","NHL","MM","HL")
dm=mkDF(pm,brkS)
df=mkDF(pf,brkS)
df=df[!df$cancer1%in%HM,]  # exclude non-heme first cancers
dm=dm[!dm$cancer1%in%HM,]
paste(levels(df$int),collapse=", ") # to paste intervals into Figure caption
# sum observed and expected cases across all non-heme first cancer types
f=df%>%group_by(trt,int)%>%summarize(O=sum(O),E=sum(E),t=mean(t,na.rm=T),sex="Female")
m=dm%>%group_by(trt,int)%>%summarize(O=sum(O),E=sum(E),t=mean(t,na.rm=T),sex="Male")
D=rbind(f,m)
(D=D%>%mutate(RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)))
D$Radiation="No" # update treatment names for legend
D$Radiation[D$trt=="rad"]="Yes"
D$Period="Early"  # make two separate time scales for plotting
D$Period[D$t>5]="Late"
graphics.off()
quartz(width=6.5,height=3.8)
theme_set(theme_bw())
theme_update(legend.position = c(.80, .35))
theme_update(axis.text=element_text(size=rel(1.4)),
             axis.title=element_text(size=rel(1.4)),
             legend.title=element_text(size=rel(1.1)),
             legend.text=element_text(size=rel(1.1)),
             strip.text = element_text(size = rel(1.4)))
g=qplot(x=t,y=RR,data=D,col=Radiation,geom=c("line","point"),
        xlab="Years Since Dx of Non-Hematologic 1st Cancer",ylab="CLL 2nd Cancer Relative Risk")
g=g+scale_y_log10(breaks=c(0.3,1,10),labels=c("0.3","1","10"),limits=c(0.2,20))
g=g+facet_grid(sex~Period,scales="free")+geom_abline(intercept=0, slope=0)
g = g + scale_color_grey(start = 0, end = 0.6)
g=g+  geom_errorbar(aes(ymin=rrL,ymax=rrU,width=.15))+ theme(legend.key.height=unit(.45, "cm"))
print(g)
ggsave("~/Results/CLL/timeCrs.eps") 
