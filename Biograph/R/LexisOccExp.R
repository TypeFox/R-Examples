LexisOccExp <-
function (Bdata,transition,nyear)
{ # NOTE: transition need origin state and destination state (hence NOT *M)
	namstates <- attr(Bdata,"param")$namstates
  if (!exists("namstates")) z<- StateSpace(Bdata)
  if (substr(transition,1,1)%in%namstates==FALSE) stop ("LexisOccExp: Origin state is not in state space")
  if (substr(transition,2,2)%in%namstates==FALSE) stop ("LexisOccExp: Destination state is not in state space")
  z <- check.par (Bdata)
  namstates <- attr(Bdata,"param")$namstates
  locpat <- locpath(Bdata)
  transition1 <- substr(transition,1,1) # State of origin (for open interval)
 # require (Epi)
	# select the subjects that experience the given transition (eg "JN") (determine their ID)
 # subjectsID1 <- Bdata$ID[grep(transition,Bdata$path,value=FALSE)]
  subjectsID1 <- TransitionAB(Bdata,transition=transition)$id
  subjects1 <- which (Bdata$ID %in%subjectsID1)
 # subjects1 <- grep(transition,Bdata$path,value=FALSE) # line number of subject with given transition
  # Select subjects that 
  #      a) experiences the transtions, OR
  #      b) do not experience the transition and are censored in the origin state (J, say)
  # Select subjects who are censored in the given origin state (J,say)
  kk <- substr(transition,1,1)
  z <- apply (Bdata,1,function (x) substr(x[locpat],x[(locpat-1)],x[(locpat-1)])==substr(transition,1,1))
  # Select subjects who do not experience given transition but are censored in origin state (determine their ID)
  subjectsID2<- Bdata$ID[(Bdata$ID %in% subjectsID1)==FALSE & z==TRUE]
  subjects2 <- Bdata$ID %in% Bdata$ID[subjects1]==FALSE & z==TRUE
  # BdataT : data for subjects with the given transition OR who are censored in origin state (J)
  BdataT <- Bdata[Bdata$ID %in%c(subjectsID1,subjectsID2),]
  ns <- nchar(BdataT$path)
  # Number of closed intervals
   print (paste("Closed intervals  =  ",length(na.omit(subjectsID1)),sep=""))
   # Number of open intervals
   print (paste("Open   intervals  =  ",length(na.omit(subjectsID2)),sep=""))

  # Determine the position of entry into risk set (first entry) [transition)))[1]]  
  # Determine the starting date and ending date of episode (exposure) in YEARS
  BdataT <- date_b(Bdata=Bdata,selectday=1,format.out="year",covs=NULL)  
    #  table(BdataT$path)
  pos <- vector (mode="numeric",length=nrow(BdataT))
  Tstart <- vector (mode="numeric",length=nrow(BdataT))
  Tstop <- vector (mode="numeric",length=nrow(BdataT))
  Tstatus <- vector (mode="numeric",length=nrow(BdataT))
  for (i in 1:nrow(BdataT))
  { if (BdataT$ID[i]%in%subjectsID1)  # closed interval
   	 { pos[i] <- nchar(unlist(strsplit(BdataT$path[i],transition)))[1]+ 1 } else
   	 { pos[i] <- nchar(unlist(strsplit(BdataT$path[i],transition1)))[1]+ 1 }
  }
  # SURVIVAL OBJECT
  
  for (i in 1:nrow(BdataT))
    { zz <- ifelse (pos[i]==1,BdataT$start[i],BdataT[i,(locpat+pos[i]-1)])
      if (BdataT$ID[i] %in% Bdata$ID[subjects1]) # closed interval
      { Tstart[i] <- zz
      	Tstop[i]  <- BdataT[i,(locpat+pos[i])]
      	Tstatus[i] <- 1 } else
      { Tstart[i] <- zz # BdataT[i,(locpat+ns[i]-1)]
      	Tstop[i] <- BdataT$end[i]
      	Tstatus[i] <- 0 }
    }
    Tstop <- ifelse (is.na(Tstop),BdataT$end,Tstop)
    
   #   Tstart[i] <- ifelse (BdataT$ID[i] %in% Bdata$ID[subjects1]== TRUE,
   #           zz,BdataT[i,(locpat+ns[i]-1)])
   #   Tstop[i] <- ifelse (BdataT$ID[i] %in% Bdata$ID[subjects1]== TRUE,
   #           BdataT[i,(locpat+pos[i])],BdataT$end[i])
   #   Tstatus[i] <- ifelse (BdataT$ID[i] %in% Bdata$ID[subjects1]== TRUE,1,0) }
  
 # ===============  Create Lexis object ======================
  print  ("Create Lexis object",quote=FALSE)
  # Transform data in years   
  bt <- BdataT$born
  endt <- BdataT$end
  en1 <-  BdataT$start
  ex1 <- Tstart
  en2 <- ex1
  ex2 <- Tstop
  event2 <- Tstatus
    duration <- ex2 - en2
  duration.neg <- length (duration[duration < 0]) # number of negative durations
  if (duration.neg > 0)
  { print ("Lexislines.episodes.R: some durations are negative.")
    print (Bdata[duration<0,])
    return
  }
  Lcoh <- Lexis( id = BdataT$ID,
               entry = list( CalTime=en2),
               exit  = list( CalTime=ex2, Age=ex2-bt ),
               exit.status = event2)

# Split Lexis object into disjoint follow-up intervals of nyear years
N <- 1   #   365.25
AgeLow <- nyear*trunc(min(na.omit(en2-bt)/nyear))
AgeHigh <- nyear*trunc(max(na.omit(ex2-bt)/nyear))
PerLow <- nyear* trunc(min(na.omit(en2)/nyear))
PerHigh <- nyear* (trunc(max(na.omit(ex2)/nyear))+1)
PerHigh[PerHigh-PerLow < AgeHigh-AgeLow] <- PerLow + AgeHigh - AgeLow
AgeHigh[AgeHigh-AgeLow < PerHigh-PerLow] <- AgeLow + PerHigh - PerLow

Lcoh_tr1_p <- splitLexis(Lcoh, breaks=seq(PerLow,PerHigh,nyear), time.scale="CalTime" )
# Lcoh_tr1_ap is Lexis object with observations on individuals divided in disjoint follow-up intervals of nyear years
Lcoh_tr1_ap <- splitLexis(Lcoh_tr1_p, breaks=seq(AgeLow,AgeHigh,nyear), time.scale="Age" )
Lcoh_tr1_ap$AGE <- timeBand(Lcoh_tr1_ap,"Age","left")
Lcoh_tr1_ap$PER <- timeBand(Lcoh_tr1_ap,"CalTime","left")
Lcoh_tr1_ap$COHORT <- (Lcoh_tr1_ap$PER - Lcoh_tr1_ap$AGE)
zzz <- Lcoh_tr1_ap$CalTime - Lcoh_tr1_ap$COHORT   #CalTime = year of birth (real value)
Lcoh_tr1_ap$UL <- ifelse (zzz < 0, 1,0) # Upper and Lower Triangle
        # EPI p. 78; APC in computer age p. 5
#       To get figures in Biographies. Real and synthetic, multiply by 365.25
# Determine number of events by Age-Period interval
nevents <- tapply (status(Lcoh_tr1_ap,"exit")==1,
  list(Lcoh_tr1_ap$AGE,Lcoh_tr1_ap$PER),sum)
#  list(timeBand(Lcoh_tr1_ap,"Age","left"),timeBand(Lcoh_tr1_ap,"CalTime","left")),sum)
# Determine exposure time by Age-Period interval   (in days)
ndur <- round(N * tapply (dur(Lcoh_tr1_ap),list(Lcoh_tr1_ap$AGE,Lcoh_tr1_ap$PER),sum),2)
ndur <- ifelse(is.na(ndur),0,ndur) # replace NA by 0 for Lexis diagram
rates <- nevents / ndur           # First transition rate PER month
date.mid <- timeBand(Lcoh_tr1_ap,"CalTime","mid")
age.mid <-  timeBand(Lcoh_tr1_ap,"Age","mid")
# =============  Draw Lexis diagram  ============================
#  EXPOSURE  (in years)
title1 <- paste ("Transition ",transition,": exposure in (years)",sep="")
print (paste("Exposure time = ",sum(ndur,na.rm=TRUE), " years or ",sum(ndur,na.rm=TRUE)*12," months",sep=""))
Lexis.diagram( age=c(AgeLow,AgeHigh), date=c(PerLow,PerHigh), coh.grid=FALSE,int=5,lab.int=5,main=title1)
date66 <- sort(unique(date.mid))
age66 <- sort(unique(age.mid))
for (ix in 1:length(age66)) {for (iy in 1:length(date66))
        text( date66[iy],age66[ix], trunc(ndur[ix,iy]), cex=0.7 ) }
par (ask=TRUE) # Display graph but hit on console
# 
# EVENTS
title1 <- paste ("Transition ",transition,": occurrences",sep="")
sum(nevents,na.rm=TRUE)
Lexis.diagram( age=c(AgeLow,AgeHigh), date=c(PerLow,PerHigh), coh.grid=FALSE,int=5,lab.int=5,main=title1)
for (ix in 1:length(age66)) {for (iy in 1:length(date66))
        text( date66[iy],age66[ix], trunc(nevents[ix,iy]), cex=0.8 ) }
par (ask=TRUE) 
#
# RATES
title1 <- paste("Transition ",transition,": Occurrence-exposure rates (per year)",sep="" )
Lexis.diagram( age=c(AgeLow,AgeHigh), date=c(PerLow,PerHigh), coh.grid=FALSE,int=5,lab.int=5,main=title1)
#
for (ix in 1:length(age66)) {for (iy in 1:length(date66))
        text( date66[iy],age66[ix], round(nevents[ix,iy]/ndur[ix,iy],3), cex=0.7 ) }
# close.screen()
par (ask=TRUE) 

# Survival object: dates in years
  surv <- Surv(Tstart,Tstop,Tstatus)
  
return (list(
			 surv=surv,
			 Lcoh = Lcoh,
             nevents=nevents,
             ndur = ndur,
             rates = rates))
}
