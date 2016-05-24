 rm(list=ls())
library (Biograph)
data (GLHS)
str(GLHS)
# Summary measures that characterize the data
param <- Parameters (GLHS)
param
# Determines the state sequences in the population under observation
seq <- Sequences(GLHS)
seq$sequences[1:5,]
# Determine state occupancies by age
occup <- Occup(GLHS)
# Display state occupancies by age
occup$state_occup
# Display, for each subject under observation, the state occupied at birthdays. This data structure is used as input in the TraMineR package.
# occup$st_age_1
# Display, for each subject under observation and age, the number of years spent in each of the states of the state space
# occup$sjt_age_1 
    # This is an array with 201 rows (ID), 54 columns (Age) and 4 layers (2 states, "censored" and "Total")
# Express the event dates in age instead of CMC
GLHS.a <- date_b (GLHS,format.out="age",covs=c("marriage","LMentry"))
# Plot the state occupancies by age
occup.a <- Occup(GLHS.a)
z<- plot.occup.S (x=occup.a$state_occup,namstates.desired=c("N","J","Censored"),colours=c("red","green","lightgrey"),title="States occupancies. GLHS",area=TRUE,xmin=10,xmax=55) 
# Get ages at transition
agetrans <- AgeTrans (GLHS)
agetrans$ages  # ages
agetrans$st_censoring  # state at censoring
table(agetrans$st_censoring)
# Censoring by age and sex
table(cut(agetrans$agecens, breaks=attr(GLHS,"param")$namage,
include.lowest=TRUE, right=FALSE), GLHS$sex)
# States at censoring by sex
table(agetrans$st_censoring,GLHS$sex)

# Tabulate number of transitions (all) by age and calendar year (p.32)
   # year of transition
   yeartrans <- YearTrans (GLHS)
a1 <- cut(agetrans$ages,breaks=seq(0,55,by=5))  # length(a1) = 2412
b1 <- cut(yeartrans[,5:ncol(yeartrans)],breaks=seq(1943,1983,by=5))  # length(b1) = 2412
z <- table (a1,b1)
z

# Tabulate transitions (all) by birth cohort (p. 33)
locpat <- locpath (GLHS)
nsample <- nrow(GLHS)
a1mat <- matrix(a1,c(nsample,ncol(GLHS)-locpat))
b1mat <- matrix(a1,c(nsample,ncol(GLHS)-locpat))

bb<- trunc(data.frame(yeartrans)$born)
c1 <- cut(bb,breaks=seq(5*trunc(min(bb)/5),5*trunc(max(bb)/5+1),by=5))
table(a1mat[,1],b1mat[,1],c1)

# Extract information on a single transition (p. 34)
z <- TransitionAB (GLHS,"NJ")
str(z)
mean(z$age,rm.na=TRUE)
round (mean(z$age,rm.na=TRUE),2)
#  Number of NJ transitions by age, cohort and calendar year
table(trunc(z$age),trunc(z$year),trunc(z$cohort))
#  Transitions by age (5yrs) and birth cohort
a1 <- cut(z$age,breaks=seq(0,55,by=5))
table (a1,c1)
# by one year: table(trunc(z$age),trunc(z$year),trunc(z$cohort))
# by 5 years
aNJ <- cut(z$age,breaks=seq(0,35,by=5)) 
yNJ <- cut(z$year,breaks=seq(1943,1983,by=5)) 
cNJ <- cut(z$cohort,breaks=seq(5*trunc(min(bb)/5),5*trunc(max(bb)/5+1),by=5))
table (age=aNJ,cohort=cNJ)
table (age=aNJ,year=yNJ,cohort=cNJ)
# Plot frequency diagram (histogram)
plot(table(trunc(z$age)),ylab="Count",xlab="Age",
   main="Number of first NJ transitions by age")

# -------  Transitions  ------------
#  Transitions by age, origin and destination
trans <- Trans (GLHS)
trans$Ttrans


# -------  ratetable ---------
z<- Parameters (GLHS)
occup <- Occup(GLHS)
trans <- Trans (GLHS)
w <- RateTable(GLHS,occup, trans)
str(w)
w$Stable

# Split GLHS by birth cohort and determines most frequent sequences for each birth cohort
GLHS.cohort <- split (GLHS,as.factor(GLHS$cohort))
str(GLHS.cohort)
# Sequences by birth cohort
seq.c1 <- Sequences (GLHS.cohort[[1]])
seq.c2 <- Sequences (GLHS.cohort[[2]])
seq.c3 <- Sequences (GLHS.cohort[[3]])
seq.c1$sequences[1:5,] 
seq.c2$sequences[1:5,]
seq.c3$sequences[1:5,]
# Display age at labour-market entry (density) by cohort and sex
library (lattice)  
z <- TransitionAB(GLHS,"NJ")
zzz <- as.data.frame(cbind (ID=GLHS$ID,
      cohort=GLHS$cohort,sex=GLHS$sex,Entry=z$age))
zzz$cohort <- factor(zzz$cohort,labels=c("1929-31","1939-41","1949-51"))
zzz$sex <- factor(zzz$sex,labels =c("Males","Females"))

# densityplot (p. 53) (lattice)
densityplot (~Entry|sex,data=zzz,plot.points="rug",
      main="Age at labour market entry",
      sub= paste("Total number of entries with known covariates is ",
      length(na.omit(zzz$Entry)),sep=""), xlab="Age", 
      scale=list(x=list(alternating=FALSE)), groups=cohort,
      ref=TRUE,auto.key=TRUE)

# Plot duration of job and no-job episodes and distinguish open and closed intervals
# First obtain an episode file (long format)
Dlong <- Biograph.long (GLHS)
str(Dlong)
# Duration of episodes (J and N) in months
Dlong2 <- Dlong$Depisode
Dlong2$time <- Dlong2$Tstop-Dlong2$Tstart
title1 <- "Duration of episodes (J and N) (months). GLHS"
subtitle1 <- "Data from cohorts 1929-31, 1939-41 and 1949-51; 982 episodes"
rownames(Dlong2) <- 1:nrow(Dlong2)
plot(rownames(Dlong2)[Dlong2$status==1 & Dlong2$OR=="J"],
     Dlong2$time[Dlong2$status==1& Dlong2$OR=="J"], 
     xlab="Episode ID",ylab="Duration of episode" , 
     pch=19, ylim=c(0,450), xlim=c(0,nrow(Dlong2)),
     main=title1,sub=subtitle1,
     cex.sub=0.6,col="red")
points(rownames(Dlong2)[Dlong2$status==0 & Dlong2$OR=="J"],        Dlong2$time[Dlong2$status==0& Dlong2$OR=="J"],pch=1,col="red")
points(rownames(Dlong2)[Dlong2$status==1 & Dlong2$OR=="N"],                    Dlong2$time[Dlong2$status==1& Dlong2$OR=="N"],pch=19,col="blue")
points(rownames(Dlong2)[Dlong2$status==0 & Dlong2$OR=="N"],Dlong2$time[Dlong2$status==0& Dlong2$OR=="N"],pch=1,col="blue")
legend (1, 460, c("Closed interval", "Open interval"),pch=c(19,1),bg="white"  )
legend (780,460, c("J episodes","N episodes"),
                    fill=c("red","blue"),cex=0.9,bg="white" )
abline (v=seq(0,1000,by=100),lty=2,col="lightgrey")
abline (h=seq(0,450,by=50),lty=2,col="lightgrey")

# Show the employment career for respondents with IDs 120 and 208
# Employment career for respondents 120 and 208
Dlong$Depisode[Dlong$Depisode$ID %in% c(120,208),]

# LEXIS diagram
# Display labour market entries by age and calendar date (scatter plot)
yeartrans <- YearTrans(GLHS)
library (Epi)
y <- Parameters(GLHS)
# Display labour market entries for total sample
z <- Lexispoints (GLHS,"NJ","Calendar time and age at labour market entry",NA)
# Display labour market entries by level of education
zE <- Lexispoints (GLHS,"NJ","Calendar time and age at labour market entry","edu")

# Linelines with on x-axis calendar year instead of CMC
   # First convert CMC to calendar years
GLHS.yr <- date_b (GLHS,format.out="year",covs=c("marriage","LMentry"))
Dlong.yr <- Biograph.long (GLHS.yr)
tit5 <- "Employment careers for a selection of subjects. GLHS"
subjects <- c(1,78,120,208)
GLHSd <- Remove.intrastate (GLHS)
z2 <- Lexislines.episodes (Bdata=GLHSd,Dlong=Dlong.yr$Depisode,subjectsID = subjects,title = tit5)

# Lexis diagram with event counts, exposure times and transition rates
  #  event = transition from job to nojob
w <- LexisOccExp (Bdata=GLHS, transition= "JN",nyear=5)

# Compute age-period rates of transitions: rates by calendar year (1940=85) and age (10-55)
PerLow <- 1940
PerHigh <- 1985
nyear <- 5
Lcoh_tr1_p <- splitLexis(w$Lcoh, breaks=seq(PerLow,PerHigh,nyear),time.scale="CalTime" )
AgeLow <- 10
AgeHigh <- 55
Lcoh_tr1_ap <- splitLexis(Lcoh_tr1_p,breaks=seq(AgeLow,AgeHigh,nyear), time.scale="Age" )
Lcoh_tr1_ap$AGE <- timeBand(Lcoh_tr1_ap,"Age","left") 
Lcoh_tr1_ap$PER <- timeBand(Lcoh_tr1_ap,"CalTime","left")
nevents <- tapply (status(Lcoh_tr1_ap,"exit")==1,list(Lcoh_tr1_ap$AGE,Lcoh_tr1_ap$PER),sum)
ndur <- round(tapply (dur(Lcoh_tr1_ap), list(Lcoh_tr1_ap$AGE,Lcoh_tr1_ap$PER),sum),2)
rates <- nevents / ndur           # First transition rate PER month
nevents
ndur 
rates

# Produce a perspective plot: event count by calendar year and age
ages <- seq(10,50,by=5)
years <- seq (1940,1980,by=5)
persp(ages, years, nevents, theta = 30,phi = 30, 
expand = 0.5, col = "green",ticktype = "detailed")

# Plot state distribution by age using package TraMineR
occup <- Occup(GLHS) 
install.packages("TraMineR")
library (TraMineR)
DTraMineR <- seqconc (occup$st_age_1,sep="-")
namstates <- attr(GLHS,"param")$namstates
namst <- c(namstates,"C")
D.seq <- seqdef (DTraMineR,states=namst)
seqplot(D.seq, type="d",title="State distribution. GLHS", ylab="Count", xtlab=0:54, group=GLHS$sex)
seqiplot(D.seq,tlim=GLHS$ID%in%c(1,20,208))
seqfplot(D.seq,group=GLHS$sex,title="Sequence frequency plot.GLHS",ylab="Percentage",xtlab=c(0:54), ltext=c("N","J","Censored"))

# ================  SURVIVAL ANALYSIS   =====================
# Survival analysis of job episodes (repeats analysis by Blossfeld and Rohwer, 2002)
Dlong <- Biograph.long (GLHS)
# Select the 600 job episodes
D <- subset (Dlong$Depisode,Dlong$Depisode$OR=="J")
# Create a survival object
library (survival)
# Duration of job episode in months: time = Tstop - Tstart
D$time <- D$Tstop - D$Tstart # time = length of job episde
surv1 <- Surv(D$time,D$status)   
# Create survival object with starting time and ending time of episodes
surv2 <- Surv(D$Tstart,D$Tstop,D$status)  

# -----------  Empirical survival function: Kaplan-Meier   ---------------
KM <- survfit (Surv(time,status)~sex,data=D,conf.type="plain")
# summary(KM)
plot(KM,conf.int=TRUE,xlab="Duration of job episode (months)",ylab="survival probability",col=c("red","blue"),
    main="Kaplan-Meier estimator of job duration, GLHS data",
    sub="Data from cohorts 1929-31, 1939-41 and 1949-51; 600 job episodes; Compare BR2002 p. 78",cex.sub=0.7)
legend(30,1,c("Males","Females"),col=c("red","blue"),fill=c("red","blue"),cex=0.9,bg="white")

# ----------   Exponential transition model   -----------------------------
# Blossfeld and Rohwer: event count/exposure time = rate: 458/40781 = 0.01123 per month 
#    (B&R, p. 86 and 93) 
D$one <- 1
tra_exp0 <- survreg(Surv(time,status)~one,data=D,dist="exponential")
tra_exp0$coefficients
sqrt(tra_exp0$var)  # error term of intercept shown by BR 2002 p. 93
print (paste("Monthly rate = ",round(exp(-tra_exp0$coefficients[1]),4),sep=""),quote=FALSE)
# R                 4.49       0.0467
# TDA  Constant    -4.4891     0.0467
# TDA    rate = exp(-4.4891 ) =0.01123
print (paste("Overall rate of leaving job =",round(exp(-tra_exp0$coefficient[1]),6),sep=""),quote=FALSE)
# Exponential is same as Weibull with scale = 1
survreg(Surv(time,status)~one,data=D,dist="weibull",scale=1)
#  predict duration in job (is equal to 1/exp(-tra_exp0$coefficients[1])  )
#             -4.4891     0.0467   in TDA
tra_exp0_pred <- predict (tra_exp0,se.fit=TRUE)
print (paste("Expected length of job episode = ",round(tra_exp0_pred$fit[1],2)," months. Standard error = ",round(tra_exp0_pred$se.fit[1],3)," months",sep=""),quote=FALSE)  
# with covariates
tra_exp_sex <- survreg(Surv(time,status)~as.factor(sex),data=D,dist="exponential")
summary(tra_exp_sex)

# 
# --------  Event count and person-years of exposure  -------------
pyears(Tstop-Tstart ~ tcut(born,c(0,1000)),data=D,scale=1)
 # results in 600 job exits and 40762 months of observation (see OVERVIEW in Surveylife)
# Event count and person-years by agezz <- pyears(Tstop-Tstart ~ cut(born,c(0,46,1000)),data=D,scale=1)  # NO
D$age <- (D$Tstart-D$born)/12
zz_age <- pyears(time ~ tcut(age,c(0,15:50,100)),data=D,scale=1)

# --------------  Cox model   --------------
one <- rep(1,600)
Cox_s <- coxph(Surv(Tstop-Tstart,status) ~ sex,data=D,method="breslow")   # reference = males
summary (Cox_s)

# Test proportional hazard assumption
Cox_s.zph <- cox.zph(Cox_s,transform="identity",global=TRUE)   # CHECK
plot (Cox_s.zph)

# Cumulative hazard: males and females
# Cumulative hazard
Cox_s <- coxph(Surv(time,status) ~ +strata(sex),data=D,method="breslow")
sfits <- survfit(Cox_s)
# z <- coxph.detail(Cox_s)
zb <- basehaz(Cox_s)
zb[c(1,106),1]
#survival curve
 plot (sfits[1],conf.int=T,lty=c(1,2,2),
 main="Probability that a job spell exceeds given durations. GLHS",
     xlab="Job duration (months)",ylab="Survival probability",col="red")
lines (sfits[2]$time,sfits[2]$surv,lty=1,col="blue")
lines (sfits[2]$time,sfits[2]$lower,lty=2,col="blue")
lines (sfits[2]$time,sfits[2]$upper,lty=2,col="blue")
legend ("topright",legend=c("Males","Females"),col=c("red","blue"),fil=c("red","blue"),
      cex=0.9,bg="white")
 plot (sfits[2],conf.int=T,lty=c(1,2,2),fun="cumhaz",
 main="Cumulative hazard function. GLHS",
     xlab="Job duration (months)",ylab="Cumulative hazard",col="blue")
lines (sfits[1]$time,-log(sfits[1]$surv),lty=1,col="red")
lines (sfits[1]$time,-log(sfits[1]$lower),lty=2,col="red")
lines (sfits[1]$time,-log(sfits[1]$upper),lty=2,col="red")
legend ("topright",legend=c("Males","Females"),col=c("red","blue"),fil=c("red","blue"),
      cex=0.9,bg="white")

# Predict job duration for person with set of attributes
D$LFX <- D$Tstart - D$LMentry
Cox_f <- coxph (Surv(Tstop-Tstart,status)~edu+as.factor(cohort)+LFX,data=D,method="breslow")
# Test proportional hazard assumption
Cox_f.zph <- cox.zph(Cox_f,transform="identity",global=TRUE)   # CHECK
plot (Cox_f.zph[1])
indiv <- data.frame(edu=8,cohort="1929-31",LFX=0)
sfit1 <- survfit(Cox_f,newdata=indiv)
plot (sfit1,las=1, xlab="Job duration (month)",ylab="Survival probability",
main="Predicted job survival for hypothetical individual based on Cox_f model. GLHS",
cex.main=0.9,conf.int=T,col="black")
zv <- cbind(colnames(indiv),t(indiv))
legend (250,1.0,legend=zv[,1],box.lty=0,cex=0.9)
legend (280,1.0,legend=zv[,2],box.lty=0,cex=0.9)

indiv <- data.frame(edu=c(9,9,19,19),cohort=c("1929-31","1939-41","1929-31","1939-41"),LFX=0)
sfit <- survfit(Cox_f,newdata=indiv)
colours <- c("black","red","blue","green")
plot (sfit,las=1,col=colours, xlab="Job duration (month)",ylab="Survival probability",
main="Predicted job survival for hypothetical individual based on Cox_full model. GLHS",
cex.main=0.9,conf.int=F)
legend (150,1.0,legend=c(
   "edu=9, cohort=1929-31,LFX=0",
   "edu=19,cohort=1939-41,LFX=0",
   "edu=9, cohort=1929-31,LFX=0",
   "edu=19,cohort=1939-41,LFX=0"),
cex=0.7,
col=colours,
fill=colours)

# =========  MULTISTATE SURVIVAL ANALYSIS: mstate   =====================
# Application: reversible Markov chain: transitions NJ and JN
# Convert data to mstate object
Dmstate <- Biograph.mstate (GLHS)
tmat <-   attr(Dmstate,"trans")
require (mstate)
events (Dmstate)
Dmstate$agem <- round ((Dmstate$marriage-Dmstate$born)/12,2)
Dmstate$Tstarta <- round((Dmstate$Tstart-Dmstate$born)/12,2)
Dmstate$Tstopa <- round((Dmstate$Tstop-Dmstate$born)/12,2)
Dmstate$Tstart <- Dmstate$Tstop <- Dmstate$pres <- Dmstate$NOJ <- Dmstate$TE <- Dmstate$marriage <- NULL

# Expand the covariates (see mstate manual)
Dcov <- expand.covs(Dmstate,c("sex","cohort","edu","agem"))
# The object Dcov contains the data for mstate modeling
# Display cumulative hazard for NJ and JN transitions
c1 <- coxph(Surv(Tstarta,Tstopa,status) ~ strata(trans),data=Dmstate,method="breslow")
fit1 <- msfit (c1,trans=tmat,vartype="aalen")
plot(fit1,las=1,xlab="Age",ylab="Cumulative transition rate",
  legend=c("NJ","JN"),legend.pos="topleft",col=c("red","blue"),
  main="Cumulative transition rates. GLHS")
# Fit Cox model wit predictors: sex, cohort and level of education
cs <- coxph(Surv(Tstarta,Tstopa,status) ~ sexFemale.1+sexFemale.2+
  + cohort1939.41.1+cohort1939.41.2+cohort1949.51.1+cohort1949.51.2
  +edu.1+edu.2 + strata(trans),data=Dcov,method="breslow")
# Predict job duration for person of female, born in 1939-41 and with 11 years of schooling
newdat <- data.frame(trans=1:2,
  sexFemale.1=c(2,0),cohort1939.41.1=c(1,0),cohort1949.51.1=c(0,0),edu.1=11,
  sexFemale.2=c(0,2),cohort1939.41.2=c(0,1),cohort1949.51.2=c(0,0),edu.2=11,
  strata=1:2)
fits <- msfit (cs,trans=tmat,vartype="aalen",newdata=newdat)
plot(fits,las=1,xlab="Age",ylab="Cumulative transition rate",
  legend=c("NJ","JN"),legend.pos="topleft",col=c("red","blue"),
  main="Cumulative transition rates. GLHS")
# Predict job duration for person of male, born in 1939-41 and with 11 years of schooling
newdat <- data.frame(trans=1:2,
  sexFemale.1=c(1,0),cohort1939.41.1=c(1,0),cohort1949.51.1=c(0,0),edu.1=11,
  sexFemale.2=c(0,1),cohort1939.41.2=c(0,1),cohort1949.51.2=c(0,0),edu.2=11,
  strata=1:2)
fits <- msfit (cs,trans=tmat,vartype="aalen",newdata=newdat)
plot(fits,las=1,xlab="Age",ylab="Cumulative transition rate",
  legend=c("NJ","JN"),legend.pos="topleft",col=c("red","blue"),
  main="Cumulative transition rates. GLHS")

# ============  MULTISTATE SURVIVAL ANALYSIS: msm  ==============
library (msm)
Dmsm <- Biograph.msm (GLHS)
transitions <- statetable.msm (state,ID,data=Dmsm)
# Estimate msm model and obtain observed and expected state probabilities 
#   Time in age (x-axis)
twoway2.q <- rbind(c(-0.025, 0.025),c(0.2,-0.2)) 
crudeinits.msm(state ~ age, ID, data=Dmsm, qmatrix=twoway2.q) 
GLHS.msm.a <- msm( state ~ age, subject=ID, data = Dmsm,
qmatrix = twoway2.q, obstype=2, control=list(trace=2,REPORT=1,abstol=0.0000005), method="BFGS")
GLHS.msm.a
# Display the matrix of transition rates (with confidence intervals)
qmatrix.msm (GLHS.msm.a)
# Get state probabilities
options (digits=3)
z<- prevalence.msm (GLHS.msm.a,times=seq(15,50,1))
str(z)
plot.prevalence.msm(GLHS.msm.a, mintime = 15, maxtime = 50,legend.pos=c(20,100))

# Estimate multistate model with age as time variable and sex and cohort as predictors
out_sc.msm <- msm( state ~ age, subject=ID, data = Dmsm, qmatrix = twoway2.q, method="BFGS", use.deriv=TRUE, exacttimes=TRUE,covariates = ~ as.factor(sex)+as.factor(cohort), control = list (trace = 2, REPORT = 1 ) )
# Get matrix of transition rates (intensities) for males and females born in 1939-41
qmatrix.msm (out_sc.msm,covariates=list("Male","1929-31"))
qmatrix.msm (out_sc.msm,covariates=list("Female","1929-31"))
hazard.msm (out_sc.msm)   # Jackson 2011 p. 9

# ==============   MULTISTATE LIFE TABLE   ====================
# Remove diagonal elements (optional)
GLHSd <- Remove.intrastate(GLHS)
Dd <- GLHSd
# Show the transitions
Trans(Dd)$Ttrans
cr <- Cumrates (irate=3,Bdata=Dd)
str(cr)
occup <- Occup(Dd)
ist <- Sequences.ind (Dd$path,namstates)
trans <- Trans (Dd)
ratetable <- RateTable(Dd,occup, trans)
ratetable$Stable[20:30,,]
rates <- Rates.ac(Stable=ratetable$Stable)
S <- MSLT.S(cr$NeAa[,,,1])
radix <- c(10000,0)
mslt <- MSLT.e (S,radix)
Sp <- plot (x=S$S, e0=mslt$e0,title=NULL,area=TRUE)




