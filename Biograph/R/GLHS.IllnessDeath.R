GLHS.IllnessDeath <-
function (GLHS)
{
# ===  CREATE DATA FRAME (wide format) WITH FIRST JOB, SECOND JOB AND NOJOB =========
#                (Illness-death model)
# Transition from first job to second job
  z<- check.par (GLHS)
  if (attr(GLHS,"format.date")!="CMC") stop("GLHS: dates should be in CMC, as in original data.")
  locpat <- locpath(GLHS)
  nsample <- nrow(GLHS)
  entry <- GLHS[,(locpat+1)]  # CMC entry labour market
  J1Jt <- GLHS[,(locpat+2)]-entry # duration of first job
  J1Jt <- ifelse (is.na(J1Jt),GLHS$end-entry,J1Jt)  # time at end of episode: new job or censoring
  J1Js <-ifelse (J1Jt==GLHS$end-entry|substr(GLHS$path,3,3)=="N",0,1) # status at end if episode: 0 = censored; 1 = event
# Transition from first job  to N: NJN cases
# Select subjects with NJN
zNJN<- rep(0,201)
zNJN[substr(GLHS$path,1,3)=="NJN"] <- 1 # N is always state before first job
JNt1 <- ifelse (zNJN==1,GLHS[,(locpat+2)]-entry,GLHS[,4]-entry)
JNs1 <- zNJN
# Determine subjects with JJN, JJJN among those with NJJ : N after JJ
zNJJN<- rep(0,nsample)
zNJJN[grep("N",substr(GLHS$path,4,20))] <- 1
zNJJN[J1Js!=1]<-0
JNt <- rep(0,nsample)
JNs <- rep(0,nsample)
  #JNt <- vector (mode="numeric",length=nrow(GLHS))
  #JNs <- vector (mode="numeric",length=nrow(GLHS))
  pN <- vector (mode="numeric",length=nrow(GLHS))
    for (i in 1:nrow(GLHS)) # nrow(GLHS))
  { if (zNJJN[i]==0) {JNt[i] <- JNt1[i] ; JNs[i] <-0} else
    {z <- stringf(GLHS[i,locpat])
     # determine position of N in any but first state
     pN[i] <- grep("N",z[3:length(z)])[1]+1
     JNt[i] <- unlist(ifelse (length(pN[i])==0,GLHS[i,4]-entry[i],GLHS[i,(locpat+pN[i])]-entry[i]))
     JNs[i] <- ifelse (length(pN[i])==0,0,1)
     }
  }
  JNs <- ifelse (JNs==1,1,ifelse(JNs1==1,1,0))
# create data frame and matrix (trans) with event numbers
age_entry <- (GLHS$LMentry-GLHS$born)/12
GLHS.tg <- data.frame(ID=GLHS$ID,J1Jt,J1Js,JNt,JNs,sex=GLHS$sex,cohort=GLHS$cohort,
     born=GLHS$born,age_entry=round(age_entry,2),edu=GLHS$edu)
  
  attr(GLHS.tg,"format.date") <- "CMC"
# get the transition matrix
  tmat <- GLHS.trans()
  attr(GLHS.tg,"param") <- attr(GLHS,"param")
  attr(GLHS.tg,"param")$tmat <- tmat 
  attr(GLHS.tg,"param")$ntrans <- 3
  attr(GLHS.tg,"param")$trans_possible <- ifelse (is.na(tmat),"FALSE","TRUE")
  attr (GLHS.tg,"param")$transitions <- NULL
  attr(GLHS.tg,"param")$nntrans <- NULL 
    
 return (GLHS.tg) 
}
