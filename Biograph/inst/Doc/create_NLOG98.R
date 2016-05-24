# NLOG98
#  SelectVariablesOG98
library(foreign)
# The following SPSS file is produced by Matsuo and Willekens (2003)
data <- data.frame(read.spss
      ("/Users/franswillekens/Documents/S/Survey/Og98/WP_Technical/SPSS/NLOG98_F_CMC2.sav",use.value.labels=TRUE))
# Add education
name1 <- "/Users/franswillekens/Documents/DATA/EUR/NL/OG98/Data_english//ORIG_FEM.SAV"
datog <- as.data.frame(read.spss(name1,use.value.labels = FALSE))
datog$OPL_HB[datog$OPL_HB==9] <- NA  # missing values

labels(data)[[2]]
# From r_og98.for 
#     1  id,nyear,wft,lft_op,burgs_op,nmarriag,ncohabit,cmcint,cmcb,
#     1  cmcleave, leavewhy,
#     1  cmcco1,cmce1co,cmcco2,cmce2co,cmcco3,cmce3co,cmcco4,cmcco5,
#     1  cmcma1,cmce1ma,cmcma2,cmce2ma,cmcma3,cmce3ma,
#     1  whymd1,whymd2,aantlev,cmc_k1,cmc_k2,cmc_k3,cmc_k4,cmc_k5, 
#     1  intdur2,las,livarr
# TO SORT CMC
attach(data)
nsample <- nrow(data)
namstates <- c("H","A","C","M","K")
ID <- c(1:nsample)
born <- data$cmcb_op 
interview <- data$cmcint+ 1 # interview END of month => interview beginning of month
# Covariates
kerk <- datog$KERKGEZ
kerk[kerk%in%c(3,4,5,6)] <- 3
kerk[kerk%in%c(7,8,9,10)] <- 4
kerk[kerk%in%c(98,99)] <- NA
kerk <- factor(kerk,labels=c("no religion","Roman Catholic","Protestant","other"))
educ <- datog$OPL_HB
cohortbreaks <- c(0,720,10000) 
cohortnames <- c("<1960","1960+") 
cohort <- cut(born,breaks=cohortbreaks,labels=cohortnames,include.lowest=TRUE ) 
# Create a data frame of dates of transition
cmca <- matrix(NA,nrow=nsample,ncol=19)
cmca[,1] = cmcleave
cmca[,2] = cmcco1
cmca[,3] = cmce1co
cmca[,4]= cmcco2
cmca[,5] = cmce2co
cmca[,6] = cmcco3
cmca[,7] = cmce3co
cmca[,8] = cmcco4
cmca[,9] = cmcco5
cmca[,10] = cmcma1
cmca[,11] = cmce1ma
cmca[,12] = cmcma2
cmca[,13] = cmce2ma
cmca[,14] = cmcma3
cmca[,15 ]= cmce3ma
cmca[,3] <- apply(cmca,1,function(x) ifelse (x[3]%in%x[-3],NA,x[3]))
cmca[,5] <- apply(cmca,1,function(x) ifelse (x[5]%in%x[-5],NA,x[5]))
cmca[,7] <- apply(cmca,1,function(x) ifelse (x[7]%in%x[-7],NA,x[7]))
cmca[,11] <- apply(cmca,1,function(x) ifelse (x[11]%in%x[-11],NA,x[11]))
cmca[,13] <- apply(cmca,1,function(x) ifelse (x[13]%in%x[-13],NA,x[13]))
cmca[,15] <- apply(cmca,1,function(x) ifelse (x[15]%in%x[-15],NA,x[15]))
cmca[,2] <- ifelse (!is.na(cmca[,2]) & !is.na(cmca[,10]) & cmca[,2]==cmca[,10],NA,cmca[,2])
cmca[,16] = cmc_k1
cmca[,17] <- ifelse (leavewhy=="Cohabition"&cmca[,2]!=cmca[,1],cmca[,1],NA)
cmca[,18] <- ifelse (leavewhy=="Marriage"&cmca[,10]!=cmca[,1],cmca[,1],NA)
cmca[,19] <- ifelse (leavewhy=="Other",cmca[,1],NA)
cmca[,2] <- ifelse (!is.na(cmca[,2]) & !is.na(cmca[,18]) & cmca[,2]==cmca[,18],NA,cmca[,2])

# Leaving parental home
#     leavehwy = 1 cohabitation
#                2 marriage
#                3 other
#                4 lives at home (censored)
table(data$leavewhy)
colnames(cmca) <- c("A","C","A","C","A","C","A","C","C","M","A","M","A","M","A","K","C","M","A")
cmc <- cmca[,-c(1,17,18,19,20)]
# Sort and rearrange the dates and obtain state sequence (path)
f <- Sequences.ind.0(cmc,namstates,absorb="K")
dates <- data.frame (f$d)
colnames(dates) <- paste ("Tr",1:ncol(dates),sep="")
# State sequence (path)
path <- as.character(f$path)
# Create Biograph object
NLOG98  <- data.frame (ID=ID,born=born,start=born,end=interview,kerk=kerk,educ=educ,cohort=cohort,YearInt=trunc(1900+interview/12),path=as.character(path),dates[,1:(max(nchar(as.character(path)))-1)],stringsAsFactors=FALSE)

attr(NLOG98,"format.date") <- "CMC"
require (Biograph)
param <- Parameters (NLOG98)
attr(NLOG98,"trans") <- param$tmat
test.Biograph <- "/Users/franswillekens/Documents/R/0 0 MAC/Package/TEST.Biograph"
save (NLOG98,file=paste (test.Biograph,"/OG98/NLOG98.RData",sep=""))

load (paste (test.Biograph,"/OG98/NLOG98.RData",sep=""))
nID <- sample (NLOG98$ID,500,replace=FALSE,prob=NULL)
nnn <- NLOG98[NLOG98$ID%in%nID,]
NLOG98 <- nnn
attr(NLOG98,"format.date") <- "CMC"
require (Biograph)
param <- Parameters (NLOG98)
attr (NLOG98,"trans") <- NULL
attr(NLOG98,"param") <- param
zz9  <- "/Users/franswillekens/Documents/R/0 0 MAC/Package/Biograph.TEST/OG98/"
setwd(zz9)
save (NLOG98,file="NLOG98.5450.RData")
og <- sample (NLOG98$ID,496,replace=FALSE)
og <- unique (c(og,c(12,15,19,5442)))
og2 <- NLOG98[NLOG98$ID%in%og,]
NLOG98 <- og2
attr(NLOG98,"format.date") <- "CMC"
attr(NLOG98,"format.born") <- "CMC"
attr(NLOG98,"param") <- Parameters (NLOG98)
save (NLOG98,file=paste (zz9,"NLOG98.RData",sep=""))


