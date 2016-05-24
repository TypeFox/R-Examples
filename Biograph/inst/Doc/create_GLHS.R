
# ==============================================================
# ==== Reads rrdat1.txt from the TDA website  ==================
# ====       Create Biograph object      =======================
# ==============================================================

rm(list=ls())
library (Biograph)
data(GLHS)
#GLHS[GLHS$ns==max(GLHS$ns),]
locpat <- locpath (GLHS)
which(GLHS[1,]==GLHS[1,"path"],arr.ind=TRUE)[2]


# Read the data
url.tda <- "http://oldsite.soziologie-blossfeld.de/eha/tda/cf_files/Data/RRDAT.1"
rrdat.1 <- as.matrix (read.table(file=url.tda),header=FALSE)
colnames(rrdat.1) <- c("ID","NOJ","TS","TF","SEX","TI","TB","TE","TM","PRES","PRES1","EDU")
rownames(rrdat.1) <-c(1:nrow(rrdat.1))
rrdat <- data.frame(rrdat.1)
# Define state space
namstates <- c("N","J")
# Reshape data frame from long format to wide format
d=reshape (rrdat,idvar="ID",timevar="NOJ",v.names=c("TS"),drop=c("PRES","PRES1"),direction="wide")

# Sample size
nsample <- nrow(d)
# Subject identification numbers
ID <- d$ID
# Dates of birth
born <-  d$TB
# Observation window
start <- born   
end <-d$TI
end <- end + rep(1,nsample) 
  # In BR censoring at end of month; in Biograph censoring at beginning of month
# Covariates
sex <- d$SEX
sex <- factor(sex,labels=c("Male","Female"))
EDU <- d$EDU
marriage <- d$TM
LMentry <- d$TE 
# birth cohort
cohortbreaks <- c(0,468,504,3000)
cohortnames <- c("1929-31","1939-41","1949-51") 
cohort <- cut(born,breaks=cohortbreaks,labels=cohortnames,include.lowest=TRUE )
namcov <- c("sex","edu","marriage","LMentry","cohort")

maxns <- 15
namcmc <- paste("cmc",1:maxns,sep="")
cmc <- matrix(NA,nrow=nsample,ncol=maxns,dimnames=list(c(1:nsample),namcmc))
ns <- vector(mode="numeric",length=nsample)  # number of states occupied before censoring = events - 1
path <- vector(mode="character",length=nsample) 
jj <- 1
path[jj] <- "NJ"
cmc[jj,1] <- rrdat[1,3]  
ns[1] <- 1
if (rrdat[1,4] < rrdat[1,6])   # TF < TI
   { ns[1] <- ns[1] + 1 
     cmc[jj,ns[1]] <- rrdat[1,4] + 1 
     path[jj] <- paste (path[jj],"J",sep="")
   } 
for (i in 2:nrow(rrdat))  
 {  if (rrdat[i,3] == rrdat[i,8])    # first job
     {jj <- jj + 1
      ns[jj] <- 1
      cmc[jj,ns[jj]] <- rrdat[i,3]
      path[jj] <- "NJ"
      if (rrdat[i,4] < rrdat[i,6]) 
           { ns[jj] <- ns[jj] + 1 
             cmc[jj,ns[jj]] <- rrdat[i,4] + 1 
           } 
  } else {if (rrdat[i,3] - rrdat[i-1,4] > 1)   # job episode at least one month
       {path[jj] <- paste (path[jj],"N",sep="")
        ns[jj] <- ns[jj] + 1
        cmc[jj,ns[jj]] <- rrdat[i,3]
        path[jj] <- paste (path[jj],"J",sep="")
       }  else    path[jj] <- paste (path[jj],"J",sep="")
  if (rrdat[i,4] < rrdat[i,6])    # TF < TI
     {ns[jj] <- ns[jj] + 1 
      cmc[jj,ns[jj]] <- rrdat[i,4] + 1 
     }   
       }
 }
nsmax <- max(ns)         # max number of events
for (jj in 1:nsample) 
    { ns[jj] <- ns[jj] + 1 # number of states occupied
      if (nchar(path[jj]) < ns[jj]) path[jj] <- paste (path[jj],"N",sep="")
    }

GLHS <- data.frame (ID,born,start,end,sex,EDU,marriage,LMentry,cohort,path,cmc[,1:nsmax])  
GLHS$path <- as.character(GLHS$path)
colnames(GLHS) <- c("ID","born","start","end",namcov,"path",paste("Tr",1:nsmax,sep=""))
attr(GLHS,"format.date") <- "CMC"
attr(GLHS,"format.born") <- "CMC"
param <- Parameters (GLHS)
attr(GLHS,"param") <- param
zz9 <- "/Users/franswillekens/Documents/R/0 0 MAC/Package/Biograph.TEST/GLHS"
setwd(zz9)
test.Biograph <- getwd()
save (GLHS,file=paste (test.Biograph,"/GLHS.RData",sep="")) 

