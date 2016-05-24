# APL.R  (makes Figure 3)
rm(list=ls()) 
library(SEERaBomb) 
library(reshape2)
library(dplyr)  
library(grid)
# the following were made earlier using SEERaBomb's mkSEER
load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
load("~/data/SEER/mrgd/popsae.RData") # loads in popsae (extended to ages 85-99)
# trim down columns to bare neccesities before burdening seerSet with more than it needs 
canc=canc%>%select(-reg,-COD,-radiatn,-histo3,-ICD9)
popsa=popsae%>%group_by(db,race,sex,age,year)%>%summarize(py=sum(py)) # sum on regs

secs=c("AML","MDS","APL","AMLti") # second cancers of interest
brks=c(0,0.25,1,2,3,4,5,6,8,10,12)
(brkS=paste0("b",paste(brks,collapse="_")))

if (0) {
  pm=seerSet(canc,popsa,Sex="male",ageStart=0,ageEnd=100) #pooled (races) male seerSet
  (pf=seerSet(canc,popsa,Sex="female",ageStart=0,ageEnd=100)) #pooled (races) female seerSet
  pm=mk2D(pm,secondS=secs) 
  (pf=mk2D(pf,secondS=secs)) # list object pf goes in and also comes out, with more on it
  
  pm=tsd(pm,brks=brks,trts=c("rad","noRad")) 
  (pf=tsd(pf,brks=brks,trts=c("rad","noRad")) )
  system.time(save(pm,file="~/Results/amlMDS/pmAPL.RData")) #~10 seconds 
  system.time(save(pf,file="~/Results/amlMDS/pfAPL.RData")) 
} else {
  load("~/Results/amlMDS/pmAPL.RData") 
  load("~/Results/amlMDS/pfAPL.RData") 
}

head(popsa)
head(canc)

m=canc%>%filter(cancer%in%secs)%>%mutate(age=agedx+0.5)%>%
  group_by(cancer,sex,age)%>%summarise(cases=n())
pops=popsa%>%group_by(sex,age)%>%summarise(py=sum(py))
head(m)
s=data.frame(sex=sort(unique(m$sex)))
c=data.frame(cancer=secs)
cs=merge(c,s)
pL=left_join(cs,pops)
d=left_join(pL,m)
d[is.na(d$cases),"cases"]=0 #join left missings where zero's should be, so fix this
d$cancer=factor(d$cancer,levels=secs)
d=d%>%mutate(py=py/1e5,incid=cases/py)
d$sex=gsub("^m","M",d$sex)
d$sex=gsub("^f","F",d$sex)
d$sex=factor(d$sex)
library(ggplot2)
theme_set(theme_gray(base_size = 10))
# theme_set(theme_bw(base_size = 10))

theme_update(legend.position = c(.18, .75),
             axis.text=element_text(size=rel(2)),
             axis.title=element_text(size=rel(2)),
             legend.title = element_blank(),
#              legend.title=element_text(size=rel(1.5)),
             legend.text=element_text(size=rel(1.5)))
graphics.off()
quartz(height=5,width=7) 
g <- ggplot(d,aes(x=age,y=incid,shape=sex,col=cancer))+geom_point(size=3.3)+ 
  labs(x="Age (years)",y=expression(paste("Cases per ",10^5," Person-Years")))+    
  scale_y_log10(limits=c(.02,100)) 
# g = g + scale_color_grey(start = 0.8, end = 0)
g
ggsave("~/Results/amlMDS/APLage.png")  
ggsave("~/Results/amlMDS/APLage.eps")  

HM=c("AML","MDS","CMML","CML","MPN","ALL","CLL","HCL","OL","NHL","MM","hodgkin")
dm=mkDF(pm,brks)
df=mkDF(pf,brks)
sapply(df,class)
d=rbind(cbind(df,Sex="Female"),cbind(dm,Sex="Male"))
d=d%>%filter(!cancer1%in%HM)%>%group_by(Sex,trt,cancer2,int)%>%summarize(O=sum(O),E=sum(E),t=weighted.mean(t,py,na.rm=T))
d=d%>%mutate(RR=O/E, rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E))
head(d)

# this is to get the seconds ordered properly and colors back like they were before using mkD instead of mkDF
d$cancer2=as.character(d$cancer2)
d$cancer2=factor(d$cancer2,levels=c("AML","MDS","APL","AMLti")) 

graphics.off()
quartz(width=6,height=3.8)
xlabNR="Years Since Dx of First Cancer Not Treated With Radiation"
xlabR="Years Since Dx of First Cancer Treated With Radiation"
# cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #black
theme_update(legend.position = c(.92, .83))
theme_update(axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.text=element_text(size=rel(1)),
             strip.text = element_text(size = rel(1.5)))
# NR=TRUE # for Figure C
# NR=FALSE  # for Figure B
for (NR in c(FALSE,TRUE)) {
  if (NR) D=d%>%filter(trt=="noRad") else D=d%>%filter(trt=="rad")
  D$t=D$t+(as.numeric(D$cancer2)-1)*0.05 # for CI visibility
  g=qplot(x=t,y=RR,data=D,col=cancer2,geom=c("line","point"),#xlim=c(-.1,24),
          xlab=ifelse(NR,xlabNR,xlabR),ylab="Relative Risk")
  g=g+facet_grid(Sex~.,scales="free")+geom_abline(intercept=1, slope=0)
#   g1 <- guide_legend("Second Cancer")
#   g=g + guides(color=g1) 
#   g=g + scale_colour_manual(values=cbPalette)
  # g = g + scale_color_grey(start = 0.8, end = 0)
  g=g+  geom_errorbar(aes(ymin=rrL,ymax=rrU,width=.15))+ theme(legend.key.height=unit(.45, "cm"))
  print(g)
  if (NR) ggsave("~/Results/amlMDS/APLnr.eps")  else
    ggsave("~/Results/amlMDS/APLr.eps") 
  if (NR) ggsave("~/Results/amlMDS/APLnr.png")  else
    ggsave("~/Results/amlMDS/APLr.png") 
}
########### Note: RR ratio  CI are too wide and unstable to plot.
