# CLLcomputing.R
rm(list=ls()) 
library(dplyr)  
library(SEERaBomb) 
library(scales)
cf=function (x) comma_format()(x)
cfc=function (x) paste0(comma_format()(x),collapse=", ")
# the following were made earlier using SEERaBomb's mkSEER
load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
load("~/data/SEER/mrgd/popsae.RData") # loads in popsa
# seerStats(canc,popsae)
head(canc)
# trim down columns to bare neccesities needed for this paper. 
canc=canc%>%select(-reg,-COD,-histo3,-ICD9,-db)
# answer reviewer's question regarding unknown survival and radiation types
svUnk=canc$surv==9999 # surv unknown = 9999 months
radUnk=canc$radiatn%in%c(8,9) # 8 and 9 are unknown
t1=table(svUnk,radUnk)
(t2=addmargins(t1))
78917/101584 # 77.7%
prop.table(t1) # over 97% both known => <3% either is unknown
# merge SLL in with CLL 
canc$cancer=as.character(canc$cancer)
canc$cancer[canc$cancer=="SLL"] ="CLL" 
canc$cancer=factor(canc$cancer)

#get CLLs in background incidence calcs
cf(canc%>%filter(cancer%in%c("CLL"))%>%group_by(cancer)%>%summarize(cases=n()) )
canc%>%filter(agedx>99,cancer%in%c("CLL"))%>%summarize(cases=n()) 
canc%>%filter(agedx<100,cancer%in%c("CLL"))%>%summarize(cases=n()) 
cf(sum(popsae$py)) # total PY at risk

#  get PY after firsts (note: we will subset on nonheme firsts after getting PY)
pm=seerSet(canc,popsae,Sex="male",ageStart=0,ageEnd=100) #pooled (races) male seerSet
pf=seerSet(canc,popsae,Sex="female",ageStart=0,ageEnd=100) #pooled (races) female seerSet
pm=mk2D(pm,secondS=c("CLL")) 
pf=mk2D(pf,secondS=c("CLL")) 
# plot2D(pm,col="gray") #  the pngs go to "~/Results/plots/pMs0e100/CLL.png"
# plot2D(pf,col="gray") #  the pngs go to "~/Results/plots/pFs0e100/CLL.png"
# to just look at the surfaces, override the default of write, e.g. use plot2D(pm,write=F)

########## START  get remaining numbers in Figure 1. 
# For Fig 1 we want all times t>0, so we set the brks vector to brks = 0
pm=tsd(pm,brks=0,trts=c("rad","noRad"),PYLong=T)  
pf=tsd(pf,brks=0,trts=c("rad","noRad"),PYLong=T) # get individual PY to make medians 
PYm=pm$L$b0$rad$PYL$'(0,100]' 
PYmN=pm$L$b0$noRad$PYL$'(0,100]' 
head(PYm) # cancer1 is the first cancer type, cancer2 the second cancer type
# year is the starting year of a PY strip of length py between ages ageL and ageR
PYf=pf$L$b0$rad$PYL$'(0,100]' #sexes are separated by default, so we have to merge them
PYfN=pf$L$b0$noRad$PYL$'(0,100]' #sexes are separated by default, so we have to merge them
PY=rbind(PYm,PYf)                                # with rbinds to get totals
HM=c("AML","AMLti","APL","MDS","CMML","CML","MPN","ALL","CLL","HCL","OL","NHL","MM","HL")
PYnh=PY%>%filter(!cancer1%in%HM)
cf(sum(PYnh$py))
cf(length(PYnh$py))
table(PYnh$cancer1)

PYN=rbind(PYmN,PYfN)
PYNnh=PYN%>%filter(!cancer1%in%HM)
cf(sum(PYNnh$py))
cf(length(PYNnh$py))

Om=pm$L$b0$rad$Obs$'(0,100]'
Of=pf$L$b0$rad$Obs$'(0,100]'
O=rbind(Om,Of)
cf(sum(O[!rownames(O)%in%HM,]))

Om=pm$L$b0$noRad$Obs$'(0,100]'
Of=pf$L$b0$noRad$Obs$'(0,100]'
O=rbind(Om,Of)
cf(sum(O[!rownames(O)%in%HM,]))

Em=pm$L$b0$rad$Exp$'(0,100]'
Ef=pf$L$b0$rad$Exp$'(0,100]'
E=rbind(Em,Ef)
cf(sum(E[!rownames(E)%in%HM,]))

Em=pm$L$b0$noRad$Exp$'(0,100]'
Ef=pf$L$b0$noRad$Exp$'(0,100]'
E=rbind(Em,Ef)
cf(sum(E[!rownames(E)%in%HM,]))

# zap out memory hog
pf$L$b0$rad$PYL$'(0,100]'=NULL
pf$L$b0$noRad$PYL$'(0,100]'=NULL
pm$L$b0$rad$PYL$'(0,100]'=NULL
pm$L$b0$noRad$PYL$'(0,100]'=NULL

########  now make CLL 2nd cancer RR time courses
brks=c(0,.1,.2,.3,0.6,1,1.5,2,2.5,3,4,5,7,10,13,16,20) 
pm=tsd(pm,brks=brks,trts=c("rad","noRad")) 
pf=tsd(pf,brks=brks,trts=c("rad","noRad")) 
system.time(save(pm,file="~/Results/CLL/pm.RData")) #~10 seconds
system.time(save(pf,file="~/Results/CLL/pf.RData")) #~10 seconds
