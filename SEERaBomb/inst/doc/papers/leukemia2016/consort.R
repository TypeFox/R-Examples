# consort.R  (figure 1) 
rm(list=ls()) 
library(scales);  #comma_format comes from here   
cf=function (x) comma_format()(x)
library(SEERaBomb) 
library(dplyr)  
load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
canc=canc%>%select(-reg,-COD,-radiatn,-histo3,-ICD9)
table(canc$seqnum)
canc=canc%>%filter(cancer!="benign")
table(canc$seqnum) #60-88 has been removed
#upgrade AML definition to include APL
# mapCancs# shows cancer definitions (i.e. that APL overwrites AML)
canc$cancer[canc$cancer=="APL"] ="AML" # overwrite back to AML
canc$cancer[canc$cancer=="AMLti"] ="AML" # overwrite other translocations and inversion AMLs back to AML
canc$cancer=factor(canc$cancer)  # removes APL and AMLti levels, as well as benign from above
levels(canc$cancer)

# first get the top three rows of numbers
cf(canc%>%group_by(db)%>%summarize(cases=n()) )
canc%>%filter(cancer%in%c("AML","MDS"))%>%group_by(cancer,db)%>%summarize(cases=n()) 

# now get numbers for background cases
aml=canc%>%filter(cancer=="AML")
aml%>%filter(agedx>99)%>%summarize(cases=n()) # 31 AMLs were over 99 at diagnosis will be left out
aml%>%filter(agedx<100)%>%summarize(cases=n()) # 69079 will be used to form background incidence

mds=canc%>%filter(cancer=="MDS")
table(mds$yrdx) # cases before 2001 are likely typos
mds=mds%>%filter(yrdx>2000) # so filter them out. 22 cases
mds%>%filter(agedx>99)%>%summarize(cases=n()) # 70 were over 99 at diagnosis will be left out
mds%>%filter(agedx<100)%>%summarize(cases=n()) # 49502 will be used to form background incidence

# next line not in a Figure but in Methods text for justification of spreading pop PY in 85+ out to 100
mds%>%filter(agedx>84,seqnum==2)%>%summarize(cases=n())/mds%>%filter(seqnum==2)%>%summarize(cases=n())# 21% 2nd MDS are >84

load("~/data/SEER/mrgd/popsae.RData") # loads in popsae (extended to ages 85-99)
cf(popsae%>%summarize(py=sum(py))) # total PY at risk are 1,788,450,864
cf(popsae%>%filter(year>2000)%>%summarize(py=sum(py))) # total PY at risk for MDS are 1,001,683,402
popsa=popsae%>%group_by(db,race,sex,age,year)%>%summarize(py=sum(py)) # sum on regs
cf(dim(canc)[1]) # 8,565,745 cancers, out of 8,677,429 total (benign+cancers)
# There are 5 more numbers in Figure 1 that we will get below, after we compute post first PY

load("~/Results/amlMDS/consort.RData") # has both pm and pf (made by amlMDScomputing)
PYm=pm0$L$b0$rad$PYL$'(0,100]' 
PYf=pf0$L$b0$rad$PYL$'(0,100]'
#sexes are separated by default, so we have to merge them
printPY=function(PYm,PYf) {
  D=rbind(cbind(PYm,sex="male"),cbind(PYf,sex="female"))
  print(cf(D%>%summarize(cases=n(),py=sum(py))))
  D=D%>%group_by(sex)%>%summarize(cases=n(),PY=sum(py),mean=mean(py),median=median(py))
  I=sapply(D,is.numeric)
  D[,I]=round(D[,I],2)
  cf(D)
}
printPY(PYm,PYf)
# 9,797,380 PY after firsts with radiation
# 1,906,425  cases

PYm=pm0$L$b0$noRad$PYL$'(0,100]' 
PYf=pf0$L$b0$noRad$PYL$'(0,100]'
printPY(PYm,PYf)
# 27,766,900 PY after firsts without radiation
#5,001,800 cases

Om=pm0$L$b0$rad$Obs$'(0,100]'
Of=pf0$L$b0$rad$Obs$'(0,100]'
head(Of,3) # rows are first cancers (see firstS arg of tsd), columns secondS specified by mk2D call
cf(apply(rbind(Om,Of),2,sum)) # 2,239 AML and 2,065 MDS cases observed

Om=pm0$L$b0$noRad$Obs$'(0,100]'
Of=pf0$L$b0$noRad$Obs$'(0,100]'
cf(apply(rbind(Om,Of),2,sum)) # 5,550 AML; 4,799 MDS cases observed

Em=pm0$L$b0$rad$Exp$'(0,100]'
Ef=pf0$L$b0$rad$Exp$'(0,100]'
cf(round(apply(rbind(Em,Ef),2,sum),1)) #1,225.2 AML and 1,430.0 MDS cases expected

Em=pm0$L$b0$noRad$Exp$'(0,100]'
Ef=pf0$L$b0$noRad$Exp$'(0,100]'
cf(round(apply(rbind(Em,Ef),2,sum),1)) #3,452.7 AML 3,755.5 MDS cases expected 
########  END getting numbers in consort diagram

