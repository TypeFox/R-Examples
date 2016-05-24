# timeCoursePlots.R  (for Figure 2 and S5)
# first run amlMDScomputing.R, then this
rm(list=ls())
library(dplyr)
library(reshape2)
library(ggplot2)
library(SEERaBomb)
graphics.off()
system.time(load("~/Results/amlMDS/pm.RData")) 
system.time(load("~/Results/amlMDS/pf.RData")) 
brks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12)
# heme malignancies (HM) excluded are
HM=c("AML","MDS","CMML","CML","MPN","ALL","CLL","SLL","HCL","OL","NHL","MM","hodgkin")
dm=mkDF(pm,brks)
df=mkDF(pf,brks)
sapply(df,class)
d=rbind(cbind(df,Sex="Female"),cbind(dm,Sex="Male"))
d=d%>%filter(!cancer1%in%HM)%>%group_by(Sex,trt,cancer2,int)%>%summarize(O=sum(O),E=sum(E),t=weighted.mean(t,py,na.rm=T))
d=d%>%mutate(RR=O/E, rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E))
head(d)


############# First Fig 2, then Fig S5, since calcs for Fig 2 are also used for Fig S5
quartz(width=7,height=4)
# theme_set(theme_bw())
theme_update(legend.position = c(.92, .825),
             axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             axis.title.x=element_text(size=rel(1.0)),
             legend.title=element_text(size=rel(1)),
             legend.text=element_text(size=rel(1)),
             strip.text = element_text(size = rel(1.5)))
xlabNR="Years Since Dx of First Cancer Not Treated With Radiation"
xlabR="Years Since Dx of First Cancer Treated With Radiation"
# for Figure 2 (NR (no radiation) = FALSE => radiation)
for (NR in c(TRUE,FALSE)) { 
  if (NR) D=d%>%filter(trt=="noRad") else D=d%>%filter(trt=="rad")
  D[D$cancer2=="MDS","t"]=D[D$cancer2=="MDS","t"]+0.05 # shift for CI visibility
  g=qplot(x=t,y=RR,data=D,col=cancer2,geom=c("line","point"),#ylim=c(.45,2.9),
          xlab=ifelse(NR,xlabNR,xlabR),ylab="Relative Risk")
  g=g+facet_grid(Sex~.,scales="free")+geom_abline(intercept=1, slope=0)
  # g = g + scale_color_grey(start = 0, end = 0.6)
  g1 <- guide_legend("Second\nCancer")
  g=g + guides(color=g1) 
  g=g+  geom_errorbar(aes(ymin=rrL,ymax=rrU,width=.15))+scale_y_continuous(breaks=1:6)
  print(g)
  if (NR) ggsave("~/Results/amlMDS/fig2ANR.eps")  else
    ggsave("~/Results/amlMDS/fig2A.eps")  
#   if (NR) ggsave("~/Results/amlMDS/fig2ANR.png")  else
#     ggsave("~/Results/amlMDS/fig2A.png")  
} # loop on NR

############ START Figure S5
system.time(load("~/Results/amlMDS/pm9.RData")) # 4 secs to load. 
system.time(load("~/Results/amlMDS/pf9.RData")) # 4 secs to load.
dm9=mkDF(pm9,"b0_0.25_0.5_0.75_1_1.5_2_2.5_3")
df9=mkDF(pf9,"b0_0.25_0.5_0.75_1_1.5_2_2.5_3")
d9=rbind(cbind(df9,Sex="Female"),cbind(dm9,Sex="Male"))
d9=d9%>%filter(!cancer1%in%HM)%>%group_by(Sex,trt,cancer2,int)%>%summarize(O=sum(O),E=sum(E),t=weighted.mean(t,py,na.rm=T))
d9=d9%>%mutate(RR=O/E, rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E))
head(d9)

D=rbind(cbind(d9,DB="SEER-9"),cbind(d,DB="SEER-18"))
sapply(D,class)
D=D%>%filter(t<3,cancer2=="MDS",trt=="rad")
head(D)
quartz(width=7,height=3.3)
theme_update(legend.position = c(.65, .81),
             axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.title=element_text(size=rel(1.1)),
             legend.text=element_text(size=rel(1.1)),
             strip.text = element_text(size = rel(1.5)))
g=qplot(x=t,y=RR,data=D,col=Sex,geom=c("line","point"),xlim=c(-.1,3),ylim=c(0,5),
        xlab=xlabR,ylab="MDS Relative Risk")
g=g+facet_grid(.~DB)+geom_abline(intercept=1, slope=0)
g+  geom_errorbar(aes(ymin=rrL,ymax=rrU,width=.15))
ggsave("~/Results/amlMDS/mds9vs18.eps")  
ggsave("~/Results/amlMDS/mds9vs18.png")  
