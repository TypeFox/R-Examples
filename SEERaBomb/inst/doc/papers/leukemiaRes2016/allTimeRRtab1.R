#allTimeRRtab1.R  (gets all-time RR in Table 1)
# first run CLLcomputing.R
rm(list=ls())
library(dplyr)
library(SEERaBomb)
system.time(load("~/Results/CLL/pm.RData")) # 1 secs to load. 
system.time(load("~/Results/CLL/pf.RData")) # 1 secs to load. 
mkci=function(x) sprintf("%.2f (%.2f, %.2f)",x[1],x[2],x[3])
dm=mkDF(pm,"b0")
df=mkDF(pf,"b0")
HM=c("AML","AMLti","APL","MDS","CMML","CML","MPN","ALL","CLL","HCL","OL","NHL","MM","HL")
df=df[!df$cancer1%in%HM,]
dm=dm[!dm$cancer1%in%HM,]
head(df)
f=df%>%group_by(trt)%>%summarize(O=sum(O),E=sum(E),sex="F")
m=dm%>%group_by(trt)%>%summarize(O=sum(O),E=sum(E),sex="M")
D=rbind(f,m)
(T1=D%>%mutate(RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)))
mkci(T1[1,5:7]) #F IR
mkci(T1[2,5:7]) #F no IR
mkci(T1[3,5:7]) #M IR
mkci(T1[4,5:7]) #M no IR

#now pool sexes
D1=D%>%group_by(trt)%>%summarize(O=sum(O),E=sum(E))
(T1=D1%>%mutate(RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)))
mkci(T1[1,4:6]) #IR
mkci(T1[2,4:6]) # no IR
