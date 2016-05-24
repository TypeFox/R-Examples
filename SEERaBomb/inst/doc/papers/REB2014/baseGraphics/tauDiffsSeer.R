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
dw=subset(d,race=="White")
summary(w<-mle2(cases~dpois(lambda=exp(c+age55/(taum+s*tauf))*py),method="Nelder-Mead",
             start=list(c=-10,taum=30,tauf=10),data=dw))
db=subset(d,race=="Black")
summary(w<-mle2(cases~dpois(lambda=exp(c+age55/(taum+s*tauf))*py),method="Nelder-Mead",
                start=list(c=-10,taum=40,tauf=15),data=db))

da=subset(d,race=="Asian")
summary(w<-mle2(cases~dpois(lambda=exp(c+age55/(taum+s*tauf))*py),method="Nelder-Mead",
                start=list(c=-10,taum=50,tauf=20),data=da))

kwmci=c(coef(wm)[2],confint(wm)[2,])
kwfci=c(coef(wf)[2],confint(wf)[2,])
kbmci=c(coef(bm)[2],confint(bm)[2,])
kbfci=c(coef(bf)[2],confint(bf)[2,])
kamci=c(coef(am)[2],confint(am)[2,])
kafci=c(coef(af)[2],confint(af)[2,])
rbind(kwmci,kwfci,kbmci,kbfci,kamci,kafci)
amci=c(coef(am)[1],confint(am)[1,])
afci=c(coef(af)[1],confint(af)[1,])
rbind(kamci,kafci,amci,afci)

