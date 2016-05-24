###############  A-bomb data using dplyr
rm(list=ls())
library(SEERaBomb)
mkAbomb() # converts lsshempy.csv and lssinc07.csv in ~\data\abomb into tables heme and sol in ~\data\abomb\abomb.db
db <- src_sqlite("~/data/abomb/abomb.db") #dplyr interface replaces RSQLite interface functions
class(db)
src_tbls(db)
(hd=tbl(db, sql("SELECT * from hemeDesc"))) # nrows not known?
class(hd)
str(hd)  #S3: holds where to find and what to get, but not actual values. Some subcomponents are S4. 
options(dplyr.print_min=50L) # set this option to see more rows
hd # 37 rows are all of them.  ?? => still didn't fetch them. 
options(dplyr.print_min=7L)

(H=tbl(db, sql("SELECT * from heme")))  # trims off columns as well as rows
(h=select(H,city:calg,py,agex:year,D,ATL,ALL,AMLtot,CML,CLL,NHL,MM))   # : is a nice aspect of select()
(gh = group_by(h, city,sex)  ) # grouping is now in print out header.
(gh=summarise(gh,totATL=sum(ATL),totALL=sum(ALL),totAML=sum(AMLtot),totCML=sum(CML),totCLL=sum(CLL),totNHL=sum(NHL),totMM=sum(MM),totPY=sum(py)))
# note ATL and CML city differences, and that CLL is low in japan
summarise(gh,totAML=sum(totAML),totCML=sum(totCML),totPY=sum(totPY)) # sums over sexes and removes the final grouping by city
summarise(h,totATL=sum(ATL),totALL=sum(ALL),totAML=sum(AMLtot),totCML=sum(CML),totCLL=sum(CLL),totPY=sum(py)) # go straight to ungrouped for totals
# 4e6/50# = 80000 is about right, i.e. ~80k survivors followed for ~50 years
gh  # ??  => still didn't fetch them. 
class(gh)  # tbl_sql and tbl but not a data.frame
(h=collect(h) ) #Collect() to bring it into R. Now it knows the number of rows
class(h) # and its class becomes tbl_df instead of tbl_sql and includes data.frame as a default. 


(H=tbl(db, sql("SELECT * from heme")))  # trims off columns as well as rows
H %>% filter(agex<1,D>4)%>%select(agex:D) #infants hit with more than 4Gy
# 1 of 3 getting childhood ALL supports the idea that pregnant women should not be x-rayed

(d<-collect(H %>% select(age,D,un4gy,py,AMLtot,CML,ALL,ATL,NHL,MM) %>% filter(un4gy==1)))
d=d %>%   # need to have it in R to do cut()
  mutate(dose=cut(D,c(-2,.02,.4,10),labels=c("low","med","high"),include.lowest=TRUE)) %>%
  mutate(age=cut(age,c(seq(0,80,20),110),labels=seq(10,90,20))) 
d=d%>%
  group_by(dose,age) 
d=d%>%summarize(py=sum(py),AML=sum(AMLtot)/py,CML=sum(CML)/py,ALL=sum(ALL)/py,ATL=sum(ATL)/py,NHL=sum(NHL)/py,MM=sum(MM)/py)
d

library(reshape2)
md=melt(d,id.vars=1:2,measure.vars=4:9,value.name = "Incidence")
md$age=as.numeric(as.character(md$age)) # make ages numeric to help ggplot connect the dots
library(ggplot2)
graphics.off()
theme_set(theme_gray(base_size = 16)) 
qplot(data=md,y=Incidence,x=age,col=dose,log="y")+facet_wrap(~variable)+geom_line()

##### A-bomb solids
ts=tbl(db, sql("SELECT * from solid"))  # ts = table of solids
(d<-collect(ts %>% select(age,marD,un4gy,py,lung,breast,prost,thyroid) %>% filter(un4gy==1)))
d$marD
d=d %>%   
  mutate(dose=cut(marD,c(-200,.02,.4,10),include.lowest=TRUE,labels=c("low","med","high"))) %>%
  mutate(age=cut(age,c(seq(0,80,20),110),labels=c(seq(10,70,20),90))) %>%
  group_by(dose,age) %>%
  summarise(py=sum(py),lung=sum(lung)/py,thyroid=sum(thyroid)/py,prost=sum(prost)/py,breast=sum(breast)/py)

md=melt(d,id.vars=1:2,measure.vars=4:7,value.name = "Incidence")
md$age=as.numeric(as.character(md$age))
qplot(data=md,y=Incidence,x=age,col=dose,log="y")+facet_wrap(~variable)+geom_line()
