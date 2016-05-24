
rm(list=ls())
library (Biograph)
data(NLOG98)
OG <- NLOG98
# Summary measures that characterize the data
param<- Parameters (OG)
tmat <- attr(OG,"trans")
# Determines the state sequences in the population under observation
sequences <- Sequences (OG)
sequences$sequences[1:10,]
z1 <- Sequences (OG[OG$cohort=="<1960",])
z2 <- Sequences (OG[OG$cohort=="1960+",])
z3 <- Sequences (OG[OG$cohort=="<1960" &OG$kerk=="no religion",])
z4 <- Sequences (OG[OG$cohort=="1960+" &OG$kerk=="no religion",])

# Age profile of first marriage
OG.c1 <- subset(OG,OG$cohort=="<1960")
attr(OG.c1,"param") <- attr(OG,"param")
attr(OG.c1,"format.date") <- attr(OG,"format.date")
attr(OG.c1,"format.born") <- attr(OG,"format.date")
z.c1 <- TransitionAB(OG.c1, "*M") 
OG.c2 <- subset(OG,OG$cohort=="1960+")
attr(OG.c2,"param") <- attr(OG,"param")
attr(OG.c2,"format.date") <- attr(OG,"format.date")
attr(OG.c2,"format.born") <- attr(OG,"format.date")
z.c2 <- TransitionAB(OG.c2, "*M")

# Individual life paths
samplepaths <- SamplePath (Bdata=OG,subjectsID=12)
samplepaths[[1]]
# 
subjectsID <- c(8,96,980,1056,1496,2883)
title1 ="Living arrangements. OG98" 
print ("Converting dates to decimal years for Lexis diagram.")
OG.yr <- date_b(Bdata=OG,selectday=1,format.out="year")  
Dlong2 <- Biograph.long (OG.yr)
require (Epi)
z<- Lexislines.episodes (Bdata=OG.yr,Dlong=Dlong2$Depisode,subjectsID=subjectsID,title=title1)

# Tabulate number of transitions (all) by age and calendar year (p.32)
   # year of transition
agetrans<- AgeTrans (OG)
yeartrans <- YearTrans (OG)
a1 <- cut(agetrans$ages,breaks=seq(0,55,by=5))  # length(a1) = 2412
b1 <- cut(yeartrans[,5:ncol(yeartrans)],breaks=seq(1945,2000,by=5))  # length(b1) = 2412
z <- table (a1,b1)
z
#  ================  TraMineR   =====================
library (TraMineR)
occup <- Occup(OG)
DTraMineR <- seqconc (occup$st_age_1,sep="-")
State.Occup <- occup$st_age_1[,] 
namst <- c(attr(OG,"param")$namstates,"-") 
z<- matrix(namst[State.Occup],nrow=nrow(State.Occup), ncol=ncol(State.Occup), dimnames=dimnames(State.Occup))
DTraMineRa <- seqconc (z,sep="-")
# Produce state sequence object
og.seq <- seqdef(State.Occup, 1:ncol(State.Occup),informat='STS', 								alphabet=c(param$namstates,"+"))
namstatest <- c("H","A","C","M","K","X") 
ids <- which(OG$ID%in%subjectsID) 
seqplot(og.seq,type="i",tlim=ids, ltext=namstatest, xtlab=c(0:54),withlegend="right")

# State distribution by age (and covariate)
seqplot(og.seq, type="d", title="State distribution. NLOG98",ylab="Count", xtlab=0:54)
seqplot(og.seq, type="d", 
			group=OG$cohort,
			title="State distribution.OG, by birth cohort",
			ylab="Count",
			xtlab=0:54)

# Extract information on a single transition (p. 34)
z <- TransitionAB (OG,"MK")
str(z)
mean(z$age,rm.na=TRUE)
round (mean(z$age,rm.na=TRUE),2)
meanage <- mean(z$age,na.rm=T)
mean(z$age[OG$cohort=="<1960"],na.r=T)
coh1 <- OG$cohort[OG$ID%in%z$id]
kerk1 <- OG$kerk[OG$ID%in%z$id]
meanages <- aggregate(z$age, list(cohort=coh1,Religion=kerk1), mean,na.rm=T)
z<- TransitionAB(OG,"*K") 
meanage <- mean(z$age,na.rm=T)
z<- TransitionAB(OG[OG$cohort=="<1960" &OG$kerk=="no religion",],"MK") 
hist(z$age)
mean(z$age,na.rm=T)
z<- TransitionAB(OG[OG$cohort=="1960+" &OG$kerk=="no religion",],"MK") 
mean(z$age,na.rm=T)
unique (OG$kerk)
z<- TransitionAB(OG[OG$cohort=="1960+" &OG$kerk=="Roman Catholic",],"MK") 
mean(z$age,na.rm=T)
# Display ages at first marriage by cohort and education level
library (lattice) 
transition <- "*M" 
z <- TransitionAB(OG,transition) # ages at leaving home for marriage 
za <- z$age
zzz <- data.frame(cbind (ID=z$id,cohort=OG$cohort[OG$ID%in%z$id],educ=OG$educ[OG$ID%in%z$id],trans=za))
zzz$cohort <- factor(zzz$cohort,labels=c("Born <1960","Born >= 1960")) 
zzz$educ <- ifelse (zzz$educ>4,5,zzz$educ) # recode 
zzz$educ <- factor (zzz$educ,labels=c("Primary","Secondary lower","Secondary higher","High"))
table(zzz$educ,zzz$cohort) 
sub1 <- paste("Total number of first marriages with known covariates is ",length(na.omit(za)),sep="")
densityplot (~trans|educ,data=zzz,plot.points="rug",main="Age at first marriage",sub= sub1,
xlab="Age", scale=list(x=list(alternating=FALSE)),groups=cohort,ref=TRUE, auto.key=TRUE)

# ===========  MULTISTATE LIFE TABLE  ===========
  seq.ind <- Sequences.ind (NLOG98$path,attr(OG,"param")$namstates)
   occup <- Occup (OG)
   trans <- Trans (Bdata=OG)
   ratetable <- RateTable (OG,occup=occup,trans=trans)
rates <- Rates.ac(Stable=ratetable$Stable)
S <- MSLT.S(rates$M) #MULTISTATE SURVIVAL FUNCTION based on o/e rates
radix <- c(100000,0,0,0,0)
e <- MSLT.e(S,radix)
z<- plot (x=S$S,e$e0,title="Multistate survival function. NLOG98",area=TRUE,order=attr(OG,"param")$namstates)

# ---------------------  S  --------------------
S <- S$S[,,1] # at birth: all live at parental home
# life expectancy: e$e0

# ============  MULTISTATE SURVIVAL ANALYSIS: mstate   ====================
OG$kerk[OG$kerk=="missing data"] <- NA
OG$kerk <- factor(OG$kerk,exclude="missing data") # exclude 

Dmstate <- Biograph.mstate (OG)
trans <- attr(Dmstate,"trans")
Dmstate$Tstarta <- round((Dmstate$Tstart-
Dmstate$born)/12,2)
Dmstate$Tstopa <- round((Dmstate$Tstop-
Dmstate$born)/12,2)
library(mstate)
events (Dmstate)
Dcov <- expand.covs(Dmstate,c("cohort","kerk"))
head(Dcov)
c1 <- coxph(Surv(Tstarta,Tstopa,status) ~strata(trans),data=Dmstate,method="breslow")
# msfit: computes subject-specific or overall cumulative transition hazards for each of the possible transitions in the multi-state model.
fit1 <- msfit (c1,trans=trans,vartype="aalen")
zz <- Remove.intrastate(NLOG98)
trans.Biograph <- transitions(zz)  #Biograph function (compare with "transitions" of mstate)
# Cumulative transition rates
title1 <- "Cumulative transition rates. Nelson-Aalen."
plot(fit1$Haz$time[fit1$Haz$trans==13],fit1$Haz$Haz[fit1$Haz$trans==13],type="l",xlab="Age",ylab="Cumulative transition rate",main=title1,las=1,xlim=c(15,40))
trans.selected <- c(1,2,3,7,10,13)
names.trans <- c("HC","HA","HM","HK","CA","CM","CK","AC","AM","AK","MC","MA","MK")
jj=0
for (i in trans.selected)
  { jj=jj+1
  	lines (fit1$Haz$time[fit1$Haz$trans==i],fit1$Haz$Haz[fit1$Haz$trans==i],col=palette(rainbow(length(trans.selected)))[jj])
  }
legend("topleft",legend=names.trans[trans.selected],lty=1,col=palette(rainbow(length(trans.selected))))
  
# --------------  Multistate survival function (mstate)  -------------------
namstates <- attr(OG,"param")$namstates
numstates <- length(namstates)
prob0 <- probtrans(fit1,direction="forward",predt=15)
S.mstate <- prob0[[1]][1:6] # MULTISTATE SURVIVAL FUNCTION FROM MSTATE
colours <- c("black","blue","yellow","purple","red")
title1 <- "State probabilities, estimated from NLOG98 "
plot(S.mstate[,1],S.mstate[,2],type="l",las=1,xlim=c(15,40),xlab="Age",ylab="State probability",main=title1)
lines (S.mstate[,1],S.mstate[,3],col=colours[2])
lines (S.mstate[,1],S.mstate[,4],col=colours[3])
lines (S.mstate[,1],S.mstate[,5],col=colours[4])
lines (S.mstate[,1],S.mstate[,6],col=colours[5])
for (ik in 1:numstates)
  {lines (rownames(S),S[,ik],col=colours[ik],lty=3) }
legend (22,1,legend=c("H","C","A","M","K"),col=colours,lty=1)
legend (27,1,legend=c("mstate","o/e"),lty=c(1,3),col="black")

# ===============  COX MODEL  ================
# ----------------  covariates: religion and cohort  ---------
ck <- coxph(Surv(Tstarta,Tstopa,status) ~
+cohort1960..10+cohort1960..13
+kerkRoman.Catholic.10+kerkRoman.Catholic.13
+kerkProtestant.10+kerkProtestant.13
+kerkother.10+kerkother.13
+strata(trans),
data=Dcov,
method="breslow")

# Woman is Roman Catholic
OGd <- Remove.intrastate(OG)
newdat <- data.frame(trans=1:13,
	cohort1960..10=c(0,0,0,0,0,0,0,0,0,0,0,0,0),
	kerkRoman.Catholic.10=c(0,0,0,0,0,0,0,0,0,1,0,0,0),  
	kerkProtestant.10=c(0,0,0,0,0,0,0,0,0,0,0,0,0),
	kerkother.10=     c(0,0,0,0,0,0,0,0,0,0,0,0,0),
	cohort1960..13=c(0,0,0,0,0,0,0,0,0,0,0,0,0),
	kerkRoman.Catholic.13=c(0,0,0,0,0,0,0,0,0,1,0,0,0),
	kerkProtestant.13=c(0,0,0,0,0,0,0,0,0,0,0,0,0),
	kerkother.13=     c(0,0,0,0,0,0,0,0,0,0,0,0,0),
	strata=1:13)
attr(Dcov,"param") <- Parameters (OGd)
	msck <-msfit(ck,newdata=newdat,trans=attr(Dcov,"param")$tmat)
probck <- probtrans(msck,direction="forward",predt=20)
probck[[1]]

# ==========  CUMULATIVE HAZARD:  mvna  =================
library(mvna)
Dmvna <- Biograph.mvna (OG)
zz = attr(Dmvna$D,"param")
Dmvna$D <- subset (Dmvna$D,Dmvna$D$from!="K")
attr(Dmvna$D,"param") =zz
attr (Dmvna$D,'format.date') <- "age"

na <- mvna(data=Dmvna$D,
			state.names=c("H","A","C","M","K"),
			tra=attr(Dmvna$D,"param")$trans_possible,cens.name=Dmvna$cens)

k<- cbind(time=na$'H A'$time,na=na$'H A'$na)
z<- cbind(Age=na$time,
			RiskSet=na$n.risk[,1],
			Trans=aperm(na$n.event,c(3,1,2))[,1,],
			Cens=na$n.cens)
x<- cbind (z[690:695,],k[690:695,])


xyplot(na,tr.choice=c("H A","H C","H M",
		"H K","A C","A M","A K","C A",
		"C M","C K","M A","M C","M K"),
       	aspect=1,
strip=strip.custom(bg="white",
       	factor.levels=c("HC","HA","HM",
		"HK","CA","CM","CK","AC","AM",
		"AK","AC","MA","MK"),
       	par.strip.text=list(cex=0.9)),
       	scales=list(alternating=1),
xlab="Age in years", 
xlim=c(10,60),
	       	ylab="Nelson-Aalen esimates")




