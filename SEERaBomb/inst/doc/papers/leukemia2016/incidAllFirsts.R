# incidAllFirsts.R (Figure 2)
rm(list=ls()) 
library(SEERaBomb)
library(ggplot2)
library(dplyr)
load("~/data/SEER/mrgd/cancDef.RData")
load("~/data/SEER/mrgd/popsae.RData") 
HM=c("AML","MDS","CMML","CML","MPN","ALL","CLL","SLL","HCL","OL","NHL","MM","hodgkin")

graphics.off()
quartz(height=4.5,width=7) 
theme_set(theme_bw())
theme_update(legend.position = c(.85, .21),
             axis.text = element_text(size = rel(1.5)),
             plot.title = element_text(size = rel(1.7)),
             axis.title = element_text(size = rel(1.5)),
             legend.text = element_text(size = rel(1.5)),
             legend.title = element_text(size = rel(1.5))      )

for (HSC in c(FALSE,TRUE)) {
  if (HSC)
    m=canc%>%filter(cancer%in%HM)%>%mutate(age=agedx+0.5)%>%
    group_by(sex,age)%>%summarise(cases=n()) else
      m=canc%>%filter(!cancer%in%HM)%>%mutate(age=agedx+0.5)%>%
    group_by(sex,age)%>%summarise(cases=n()) 
  
  pops=popsae%>%group_by(sex,age)%>%summarise(py=sum(py))
  head(m)
  s=data.frame(sex=sort(unique(m$sex)))
  pL=left_join(s,pops)
  head(pL)
  d=left_join(pL,m)
  d=d%>%mutate(py=py/1e5,incid=cases/py)
  head(d)
  names(d)[1]="Sex"
  levels(d$Sex)=c("Male","Female")
  g=qplot(age,incid,col=Sex,data=d,
        main=ifelse(HSC,"Hematological Cancers","Non-Hematological Cancers"),
        ylab="Cases/100,000 Person-Years",
        xlab="Age",log="y")+geom_line() + geom_vline(xintercept = c(40,50,65),col="gray")
  g = g + scale_color_grey(start = 0, end = 0.6)
  if (HSC) ggsave("~/Results/amlMDS/allCancersCombVsAgeHSC.eps") else
    ggsave("~/Results/amlMDS/allCancersCombVsAge.eps")
  if (HSC) ggsave("~/Results/amlMDS/allCancersCombVsAgeHSC.png") else
    ggsave("~/Results/amlMDS/allCancersCombVsAge.png")
}

