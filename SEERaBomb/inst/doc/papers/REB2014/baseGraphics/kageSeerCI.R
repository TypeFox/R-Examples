# NOTE: This script assumes that the SEER data is in /data/SEER 
.seerHome="/data/SEER" 
rm(list=ls()) # note: dot variables defined above persist through such cleanings
# install.packages("RSQLite")
# install.packages("bbmle")
# install.packages("ggplot2")
library(RSQLite)
m=dbDriver("SQLite")
library(plyr)
library(bbmle)
library(ggplot2)

# SEER CML versus Race Figure 4 (show sex CI separately)
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
pops=dbGetQuery(con, "SELECT * from pops where popyear>2000 and popage>5 and popage<19")
head(pops,20)
summary(as.factor(pops$poprace))
pops$poprace[pops$poprace>2]=3
(pop<-ddply(pops, .(popage,popsex,poprace), summarise,py=sum(population)))
head(pop)
d=dbGetQuery(con, 
# ICDO3 CML codes 9875 = bcr-abl+ and 9876 =bcr-abl neg CML are not in full use yet, i.e.
# many ICD-O3 CML codes are still 9863 carried over from ICD-O2
# "SELECT * from lymyleuk where histo3=9876 and seqnum<2 and race<98 and agerec>5 and agerec<19")
"SELECT * from lymyleuk where histo2=9863 and seqnum<2 and race<98 and agerec>5 and agerec<19")
d$race[d$race>2]=3
head(d)
(ddply(d, .(agerec,sex,race), summarise,cases=length(agerec))) #only 5 of the 16s 
count(d, vars=c("agerec","sex","race"))  # same here, and no options to fix it
(d=ddply(d, .(agerec,sex,race), summarise,cases=length(agerec),.drop=F))# this gets them all
head(cbind(d,pop)) # see alignment
d=cbind(d,py=pop[,"py"])
d=transform(d,incid=1e6*cases/py)
head(d,15)
d$Sex=factor(d$sex,labels=c("Male","Female"))
d$race=factor(d$race,labels=c("White","Black","Asian"))
age=c(0.5,3,seq(7.5,87.5,5))
d$age=age[d$agerec+1]
d$age55=d$age-55
d$s=d$sex-1
head(d,15)
dwm=subset(d,race=="White"&Sex=="Male")
dwf=subset(d,race=="White"&Sex=="Female")
summary(wm<-mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
             start=list(c=-12,k=0.04),data=dwm))
summary(wf<-mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
                 start=list(c=-12,k=0.04),data=dwf))
dbm=subset(d,race=="Black"&Sex=="Male")
dbf=subset(d,race=="Black"&Sex=="Female")
summary(bm<-mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
                 start=list(c=-12,k=0.04),data=dbm))
summary(bf<-mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
                 start=list(c=-12,k=0.04),data=dbf))
dam=subset(d,race=="Asian"&Sex=="Male")
daf=subset(d,race=="Asian"&Sex=="Female")
summary(am<-mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
                 start=list(c=-12,k=0.04),data=dam))
summary(af<-mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
                 start=list(c=-12,k=0.04),data=daf))
kwmci=c(coef(wm)[2],confint(wm)[2,])
kwfci=c(coef(wf)[2],confint(wf)[2,])
kbmci=c(coef(bm)[2],confint(bm)[2,])
kbfci=c(coef(bf)[2],confint(bf)[2,])
kamci=c(coef(am)[2],confint(am)[2,])
kafci=c(coef(af)[2],confint(af)[2,])
outs=rbind(kwmci,kwfci,kbmci,kbfci,kamci,kafci)
row.names(outs)<-c("White Males","White Females","Black Males","Black Females",
                   "Asian Males","Asian Females")
(outs=round(outs,4))
library(hwriter)
setwd("/users/radivot/downloads/sachs")
hwrite(outs,"table1.html")
