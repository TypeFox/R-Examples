#amlMDScomputing.R
rm(list=ls())  
library(SEERaBomb)
library(dplyr)  #the %>% operators below come from here
# the following was made earlier using SEERaBomb's mkSEER
load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
canc=canc%>%filter(cancer!="benign")
load("~/data/SEER/mrgd/popsae.RData") # loads in popsae (extended to ages 85-99)
# trim down columns to bare neccesities needed for this paper. 
canc=canc%>%select(-reg,-COD,-radiatn,-histo3,-ICD9)
# canc=canc%>%select(-reg,-recno,-agerec,-numprims,-COD,-age19,-radiatn,-histo3,-ICD9)
popsa=popsae%>%group_by(db,race,sex,age,year)%>%summarize(py=sum(py)) # sum on regs
head(canc,1)
head(popsa,1)
levels(canc$cancer) # default is that AMLti and APL are separate from AML
#next chunk upgrades AML definition to include APL and AMLti (AML by translocation or inversion)
# mapCancs# shows ICD9/ICDO3 cancer definitions 
canc$cancer[canc$cancer=="APL"] ="AML" # overwrite back to AML
canc$cancer[canc$cancer=="AMLti"] ="AML" # overwrite back to AML
canc$cancer=factor(canc$cancer)#eliminate extra levels: this happens in seerSet() too so not critical
levels(canc$cancer)

pm=seerSet(canc,popsa,Sex="male",ageStart=0,ageEnd=100) #pooled (races) male seerSet
pf=seerSet(canc,popsa,Sex="female",ageStart=0,ageEnd=100) #pooled (races) female seerSet
pm=mk2D(pm,secondS=c("AML","MDS")) # 4 secs
pf=mk2D(pf,secondS=c("AML","MDS"))# list object pf goes in and also comes out, with more on it
# comparing objects, note that mk2D added the dataframe D (for 2D plotting) to the object 

##########  background incidence plots
# plot2D(pm)  # this uses rgl to interactively pick a nice angle for a picture in FIGURE 3
# plot2D(pf,col="gray") # note from its help page in SEERaBomb that the pngs go to "~/Results/plots".
# To just look at the surfaces, override the default of write, e.g. use plot2D(pm,write=F)
##########  END background incidence

########## first some calcs over all times for the consort diagram
brks=c(0) # we want all times t>=0, so we set the brks vector to brks = 0
pm0=tsd(pm,brks=brks,trts=c("rad","noRad"),PYLong=TRUE) 
pf0=tsd(pf,brks=brks,trts=c("rad","noRad"),PYLong=TRUE) # PYLong adds a lot of bulk
system.time(save(pm0,pf0,file="~/Results/amlMDS/consort.RData")) # so save to a separate file 
######### now get full time courses for Figure 2   # this takes some computing time
brks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12) # get 12 to line up with tables
pm=tsd(pm,brks=brks,trts=c("rad","noRad")) 
pf=tsd(pf,brks=brks,trts=c("rad","noRad"))# this adds the list L to the object pf
# The list L has only one element now, named according to brks, but this will grow out
# with each call to tsd (time since diagnosis: see SEERaBomb's tsd help page) below

################### END computing RR, O, E and PY for full time courses

########### now Fig 4 first cancers less prevalent than breast or prostate
brks=c(0,1,2,3,6,9,12)   # intermediate level resolution
pm=tsd(pm,brks=brks,trts=c("rad","noRad")) 
pf=tsd(pf,brks=brks,trts=c("rad","noRad")) 
########### End specific first cancer calcs

######### START get risks over 1 to 12 years for Tables ######## 
brks=c(0,0.25,1,12)
pm=tsd(pm,brks=brks,trts=c("rad","noRad")) 
pf=tsd(pf,brks=brks,trts=c("rad","noRad")) 
system.time(save(pm,file="~/Results/amlMDS/pm.RData")) #~10 seconds 
system.time(save(pf,file="~/Results/amlMDS/pf.RData")) 

# brks=c(0,0.25,1,12)
# pm=tsd(pm,brks=brks,trts=c("rad","noRad"),PYLong=T) # this TRUE bulks things up, slows reading in 
# pf=tsd(pf,brks=brks,trts=c("rad","noRad"),PYLong=T) 
# system.time(save(pm,file="~/Results/amlMDS/pmPYM.RData")) #so save in separate file 
# system.time(save(pf,file="~/Results/amlMDS/pfPYM.RData")) 
######## END Table calcs



######### START get SEER9 time courses for the left panel of Figure S5 
canc9=canc%>%filter(db=="73")  
popsa9=popsa%>%filter(db=="73")
(pm9=seerSet(canc9,popsa9,Sex="male",ageStart=0,ageEnd=100)) #pooled (races) male seerSet
(pf9=seerSet(canc9,popsa9,Sex="female",ageStart=0,ageEnd=100)) #pooled (races) female seerSet
pm9=mk2D(pm9,secondS=c("AML","MDS")) # 3 secs, was 170 secs with all cancerS
pf9=mk2D(pf9,secondS=c("AML","MDS"))
# plot2D(pm9,write=F)  
# plot2D(pf9,write=F)
brks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3) # this is to show the improvement in variance of MDS
pm9=tsd(pm9,brks=brks,trts=c("rad")) 
pf9=tsd(pf9,brks=brks,trts=c("rad"))  
system.time(save(pm9,file="~/Results/amlMDS/pm9.RData")) #~6 seconds 
system.time(save(pf9,file="~/Results/amlMDS/pf9.RData")) 
################### END computing RR, O, E and PY for SEER9 for right panel of Fig S5. 

# the next chunk is to validate the codes via Table 1
pm20=seerSet(canc9,popsa9,Sex="male",ageStart=20,ageEnd=85) #pooled (races) male seerSet
pf20=seerSet(canc9,popsa9,Sex="female",ageStart=20,ageEnd=85) #pooled (races) female seerSet
pm20=mk2D(pm20) # 55 secs (goes through all second cancers by default)
pf20=mk2D(pf20) #
pm20=tsd(pm20,brks=c(5)) 
pf20=tsd(pf20,brks=c(5)) 
system.time(save(pm20,file="~/Results/amlMDS/pm20.RData")) #~6 seconds 
system.time(save(pf20,file="~/Results/amlMDS/pf20.RData")) 

# APL
fisher.test(matrix(c(5,13,0,13),ncol=2))

