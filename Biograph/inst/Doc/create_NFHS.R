# Andhra Pradesh NFHS3  2006

rm(list=ls())
library (foreign)
library (Biograph)
# ------------------  Read data  -------------------
zz1 <- "/Users/frans/Documents/R/India/NFHS/NFHS-3/Women/"
zz2 <- "/Users/frans/Documents/R/0 0 MAC/Package/Biograph.TEST/Chapters/AnnexA/NFHS/"
AP3 <- paste(zz2,"AP.sav",sep="")
# zz2 <- "/Users/franswillekens/Documents/R/India/NFHS2_1/"
# ka12 <- paste(zz2,"Ka2006.sav",sep="")
d06 <- data.frame(read.spss (AP3,use.value.labels=F)) 
# v011	date of birth (all dates in CMC)
# v008	date of interview 
# v509	date of first marriage
# v201	Number of children ever born (nCEB)
# b3.*	date of birth of child (from youngest to oldest)
# bord.* birth order of child
# v312	contraceptive method (sterilization = 6 (female) or 7 (male))
# v317	date of sterilization
CEB <- d06$v201 # Number of children ever born
#zz4 <- "/Users/franswillekens/Documents/R/India/NFHS/NFHS0506/IAIR52DT/"
#library(foreign)
#d.st <- read.dta (paste(zz4,"IAIR52FL.dta",sep=""),convert.dates=TRUE,convert.underscore=TRUE)
# DATA FRAME


# -----------------  Arrange births by birth order  --------
# Birth sequence (first child first). In the raw data, the
# youngest child (last birth) is listed first. Using data on birth order, the variable cmc_k06 is created, where the oldest is listed first.
cmc_child <- array(NA,dim=c(nrow(d06),20))
colnames (cmc_child) <- c(paste("ch",1:20,sep=""))
cmc_child <- cbind(d06$b3.01,d06$b3.02,d06$b3.03,d06$b3.04,d06$b3.05,
d06$b3.06,d06$b3.07,d06$b3.08,d06$b3.09,d06$b3.10,d06$b3.11,d06$b3.12,d06$b3.13,
d06$b3.14,d06$b3.15,d06$b3.16,d06$b3.17,d06$b3.18,d06$b3.19,d06$b3.20)

# Birth order
bord_child <- array(NA,dim=c(nrow(d06),20))
colnames (bord_child) <- colnames(cmc_child)
bord_child <- cbind(d06$bord.01,d06$bord.02,d06$bord.03,d06$bord.04,d06$bord.05,
d06$bord.06,d06$bord.07,d06$bord.08,d06$bord.09,d06$bord.10,d06$bord.11,
d06$bord.12,d06$bord.13,d06$bord.14,d06$bord.15,d06$bord.16,d06$bord.17,
d06$bord.18,d06$bord.19,d06$bord.20)

cmc_k06 <- array(NA,dim=c(nrow(d06),max(CEB)))
colnames(cmc_k06) <- colnames(cmc_child)
#  cmc_k06  cmc at birth of child, from first to latest birth
for (i in 1:nrow(d06)) 
  {  for (j in 1:CEB[i]){ cmc_k06[i,bord_child[i,j]] <- cmc_child[i,j]}} 
  
# --------------------  Sterilization  --------------------
cmc_ster <- vector(mode="numeric",length=nrow(d06))
for (i in 1:nrow(d06))
 { if (d06$v312[i]==6 | d06$v312[i]==7) cmc_ster[i] <- d06$v317[i] else
    cmc_ster[i] <- NA
 }

# ===============  Create Biograph object  ============================
nsample <- nrow(d06)							# sample size
namstates <- c("H","M",letters[1:max(CEB)],"S") # state space

ID <- c(1:nsample)  							# identification number
born <- d06$v011								# date of birth
start <-born									# onset of observation
end <-d06$v008 +1 # CMC at interview			# end of observation

nn <- max(CEB)+2
# -------------  dates of transition (CMC)  -----------------
cmc <- array(NA,dim=c(nsample,nn))			
cmc[,1] <- d06$v509 # first marriage
cmc[,(2:(nn-1))] <- cmc_k06
cmc[,nn]<- cmc_ster    #  d06$v317 #cmc at sterilization
dimnames(cmc) <- list(ID=ID,Transition=c("M",letters[1:max(CEB)],"S"))
cmc <- data.frame(cmc)
# ------------  transitions in chronological order  ----------
require (Biograph)
f<- Sequences.ind.0 (cmc,namstates)
path <- as.character(f$path)
# ------------  covariates  -------------------
namcov <- c("COH","EDU","WEAL","U_R")
namcohort <- c("<1970","1970-79",">=1980")
# d06$v010 is year of birth
COH <- cut(d06$v010,breaks=c(0,1969,1979,3000),include.lowest=T,labels=namcohort )
EDU <- as.factor(d06$v106)  # education
WEAL <- as.factor(d06$v190) # wealth index
U_R <- as.factor(d06$v102)  # urban/rural residence
# ------------   Make data frame  --------------
namtrans <- paste("Tr",1:ncol(f$d),sep="")
D.AP <- data.frame (ID,born,start,end,COH,EDU,WEAL,U_R,CEB,path,f$d)  
D.AP$path <- as.character(D.AP$path)
namcov <- c("COH","EDU","WEAL","U_R","CEB")
colnames(D.AP) <- c("ID","born","start","end",namcov,"path",namtrans)
locpat <- locpath(D.AP)
AP <- cbind(D.AP[,1:locpat],round(D.AP[,(locpat+1):ncol(D.AP)],2)) 
attr(AP,"format.date") <- "CMC"
attr(AP,"format.born") <- "CMC"
param <- Parameters(AP)
attr(AP,"param") <- param
save(AP,file="AP.RData")


