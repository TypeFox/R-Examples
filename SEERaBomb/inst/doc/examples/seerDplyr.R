# using dplyr with SEER data
rm(list=ls()) 
library(dplyr)
db <- src_sqlite("~/data/SEER/mrgd/cancDef.db")
src_tbls(db)    #see what tables are in it
tbl(db, sql("SELECT * from canc")) # see field names of canc
tbl(db, sql("SELECT * from popsa")) # and popga

# If you want all the rows and columns, it's faster to just load() the data.frame binary
# system.time(d<-collect(tbl(db, sql("SELECT * from canc")))) #17 secs

# db2 <- src_sqlite("~/data/SEER/mrgd/cancAll.db")
# system.time(d<-collect(tbl(db2, sql("SELECT * from canc")))) #123 secs, 123/17= 7.2
system.time(load("~/data/SEER/mrgd/cancDef.RData")) #3 secs Here Def stands for default picks in pickFields
# system.time(load("~/data/SEER/mrgd/cancAll.RData")) #30 secs
# The latter also has the advantage of keeping factor level orders. The DBs store them as strings. 

# for smaller tables and subsets of large tables, the times are negligible
collect(tbl(db, sql("SELECT * from popsa")))
# (2011-1973)*18*3*2*19 # =77976 > 54009 rows ... conditions have no person years. py=0 => no row
(mds1<-collect(tbl(db, sql("SELECT db,reg,race,sex,agedx,yrdx,COD,surv from canc where yrdx>2000 and histo3>9979 and histo3<9990"))))
# 40k MDS cases   DB was indexed on histO3 so it is fast (0.1 secs) to find the MDS cases
# Meanwhile, while
canc  # is slow to print (because it is so big)

# getting MDS cases from it  is also  fast
mds2=canc%>%
#   mutate(age=age19)%>%
  filter(yrdx>2000,histo3>9979,histo3<9990)%>%
  select(db,reg:agedx,yrdx,COD,surv,age86)
mds2
# and it has the advantage of using my original factor level orderings
table(mds2$race)
table(mds1$race)


# 1) Plot MDS incidence versus year, sex, race, and reg

# a) make a table of interest

# first group cases into 86, mostly 1-year, age buckets
mds=mds2%>%
  mutate(year=yrdx)%>%  #let year=yrdx to match population column below
  group_by(reg,db,race,sex,age86,year) %>%                                     
  summarise(cases=n())
mds
sum(mds$cases)

# now get person years
load("~/data/SEER/mrgd/popsa.RData") 
(p=popsa%>%filter(year>2000))
(d=left_join(p,mds)) 
d[is.na(d$cases),"cases"]=0 #join left missings where zero's should be, so fix this
sum(d$cases) # 40k mds cases, check.
(d1=d%>%group_by(year)%>%summarise(cases=sum(cases),py=sum(py),Incid=cases/py))
library(ggplot2)
graphics.off() # kill any left over windows
qplot(year,Incid,data=d1,ylim=c(0,6e-5),ylab="MDS Incidence") 

d=d%>%filter(year>=2006) #use only stable portion of the data. 

# look at variation across registries
(d1=d%>%group_by(db,reg)%>%
   summarise(cases=sum(cases),py=sum(py),Incid=cases/py,
             lo=Incid-1.96*sqrt(cases)/py,hi=Incid+1.96*sqrt(cases)/py ))

theme_set(theme_gray(base_size = 12)) 
qplot(reg,Incid,data=d1,ylim=c(0,8e-5),col=db,ylab="MDS Incidence (Cases/Person-Year)",xlab="SEER Registry") + 
  geom_errorbar(aes(ymin=lo,ymax=hi)) 
library(SEERaBomb)
mapRegs()

sort(table(canc$reg),decr=T)

# look at incidence vs. age for sexes and races
(d1=d%>%group_by(sex,race,age86)%>%summarise(Incid=sum(cases)/sum(py) ))
theme_set(theme_gray(base_size = 22)) 
qplot(age86,Incid,col=sex,shape=race,data=d1,log="y",ylab="MDS Incidence (Cases/Person-Year)",xlab="Age (years)")+geom_line()
# too much confusion
# trying fewer bins below doesn't help
(d1=d%>%mutate(agec=cut(age86,c(0,50,70,80,100)))%>%group_by(sex,race,agec)%>%summarise(Incid=sum(cases)/sum(py) ))
d1$agem=c(25,60,75,85)[d1$agec]
theme_set(theme_gray(base_size = 22)) 
qplot(agem,Incid,col=sex,shape=race,data=d1,log="y",ylab="MDS Incidence (Cases/Person-Year)",xlab="Age (years)")+geom_line()

with(canc,table(sex,race))

# look at survival vs. sex, age and race
library(survival)
(mds=mds2%>%filter(surv<1000)) # 1000=83 yrs. This is to get rid of some entries of 9999=survival unknown
mds=mds%>%filter(yrdx>=2006)%>%mutate(agec=cut(age86,c(0,50,70,80,100)))
coxph(Surv(surv,COD>0)~sex,data = mds) 
coxph(Surv(surv,COD>0)~sex+agec,data = mds) 
coxph(Surv(surv,COD>0)~sex+agec+race,data = mds) 
(cm=coxph(Surv(surv,COD>0)~sex+agec+race+reg,data = mds))  #best off being treated in the bay area

