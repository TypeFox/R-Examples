##desc is descriptive statistics function
desc<-function (y, pcts = c(0.025, 0.05, 0.95, 0.975), nsig = 4) 
{
    x <- y[!is.na(y)]
    rn <- range(x)
    bin <- ifelse(is.na(y), 1, 0)
    NAs <- sum(bin)
    if (NAs == 0) {
        tmp <- as.numeric(signif(cbind(N = length(x), mean = mean(x), 
            med = median(x), s = sqrt(var(x)), t(quantile(x, 
                pcts)), min = rn[1], max = rn[2]), nsig))
        names(tmp) <- c("N", "Mean", "Med", "S", paste(pcts), 
            "Min", "Max")
    }
    else {
        tmp <- as.numeric(signif(cbind(N = length(x), NAs = NAs, 
            mean = mean(x), med = median(x), s = sqrt(var(x)), 
            t(quantile(x, pcts)), min = rn[1], max = rn[2]), 
            nsig))
        names(tmp) <- c("N", "NA", "Mean", "Med", "S", paste(pcts), 
            "Min", "Max")
    }
    tmp
}

## pava is the pool adjacent violator algorithm to perform isotonic transformation
pava <- function (x, wt = rep(1, length(x))) 
{
  n <- length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}

# Simulation of clinical trial results based on toxicity
teqrOCtox<-function(sim,firstdose=2,probt,cohortSize=3,MaxNoCohorts=30,MTDss=12,pTarget,eq1,eq2,tootoxic){
results.old <- rep(NA,6)
for (j in 1:sim){
#print(j)
tox.old<-0
cumtox.old<-0
newdose<-firstdose
stopdose<-length(probt)
#print(newdose)
doselevel.old<-0
for (i in 1:MaxNoCohorts){
tox=rbinom(1,cohortSize, prob=probt[newdose])
tox.old<-rbind(tox.old,tox)
doselevel.old<-rbind(doselevel.old,newdose)
output<-data.frame(cbind(doselevel.old,tox.old), row.names=NULL)
colnames(output)<-c('doselevel','tox')
output<-output[c(-1),]
#print(paste('i=', i))
#print(output)
currentlevel<-output[output$doselevel==newdose,]
cumtox<-sum(currentlevel$tox)/(cohortSize*length(currentlevel$tox))
cumtox.old<-rbind(cumtox.old,cumtox)
#print(paste('newdose', newdose))
#print(paste('cumtox',cumtox))
upperlimit<-pTarget+eq2
lowerlimit<-pTarget-eq1
stopdosem1<-stopdose-1
#print(paste('stopdosem1', stopdosem1))
if (newdose<stopdose & newdose>1){
if (cumtox>=lowerlimit) doselevel<-newdose
if (cumtox>upperlimit) doselevel<-newdose-1
if (cumtox>=tootoxic)  stopdose<-newdose-1
if (cumtox<lowerlimit)  doselevel<-newdose+1 
}
if (newdose==1){
if (cumtox>upperlimit) break
if (cumtox<=upperlimit)  doselevel<-newdose 
if (cumtox<lowerlimit)  doselevel<-newdose+1 
}
if (newdose==stopdose){
if (cumtox>upperlimit & newdose>1) doselevel<-newdose-1
if (cumtox>=tootoxic)  stopdose<-newdose-1
if (cumtox<=upperlimit)   doselevel<-newdose 
}
newdose<-doselevel
if (cohortSize*length(currentlevel$tox)>=MTDss & cumtox<tootoxic) break
}
cumtox.m1<-cumtox.old[-c(1)]
simNo<-rep(j, length(output$tox))
results<-cbind(simNo,output, cumtox.m1)
results.old<-rbind(results.old, results)
}
simresults<-data.frame(results.old[-c(1),], row.names=NULL)
simData<-(list(simresults=simresults, cohortSize=cohortSize,probt=probt, MTDss=MTDss, pTarget=pTarget, lowerlimit=lowerlimit,upperlimit=upperlimit,tootoxic=tootoxic, sim=sim))
DLdata<-teqr.DLdata(simData=simData)
MTDdata<-teqr.MTDdatatox(simData=simData, DLdata=DLdata)
OperChartox<-teqr.OperChartox(simData=simData, DLdata=DLdata,MTDdata=MTDdata)
OperChartox
}


# Simulation of clinical trial results based on activity constrained by toxicity
teqrOCact<-function(sim, firstdose=2, proba,probc=probc, cohortSize=3, MaxNoCohorts=30, RP2Dss=RP2Dss, pTarget, eq1, eq2, tootoxic=1.01, toxcon,llactivity){
MTDss<-RP2Dss
results.old<-rep(NA,6)
for (j in 1:sim){
#print(j)
tox.old<-0
toxC.old<-0
cumtox.old<-0
cumtoxC.old<-0
newdose<-firstdose
stopdose<-length(proba)
#print(newdose)
doselevel.old<-0
for (i in 1:MaxNoCohorts){
tox=rbinom(1,cohortSize, prob=proba[newdose])
toxC=rbinom(1,cohortSize, prob=probc[newdose])
tox.old<-rbind(tox.old,tox)
toxC.old<-rbind(toxC.old,toxC)
doselevel.old<-rbind(doselevel.old,newdose)
output<-data.frame(cbind(doselevel.old,tox.old,toxC.old), row.names=NULL)
colnames(output)<-c('doselevel','tox','toxC')
output<-output[c(-1),]
#print(paste('i=', i))
#print(output)
currentlevel<-output[output$doselevel==newdose,]
cumtox<-sum(currentlevel$tox)/(cohortSize*length(currentlevel$tox))
cumtox.old<-rbind(cumtox.old,cumtox)
cumtoxC<-sum(currentlevel$toxC)/(cohortSize*length(currentlevel$toxC))
cumtoxC.old<-rbind(cumtoxC.old,cumtoxC)
#print(paste('newdose', newdose))
#print(paste('cumtox',cumtox))
#print(paste('cumtoxC',cumtoxC))
upperlimit<-pTarget+eq2
lowerlimit<-pTarget-eq1
if (newdose<stopdose & newdose>1 & cumtoxC>=toxcon){
    doselevel<-newdose-1
    stopdose<-newdose-1
}
if (newdose<stopdose & newdose>1 & cumtoxC<toxcon){
if (cumtox>=lowerlimit) doselevel<-newdose
if (cumtox>upperlimit) doselevel<-newdose-1
if (cumtox>=tootoxic)  stopdose<-newdose-1
if (cumtox<lowerlimit)  doselevel<-newdose+1 
}
if (newdose==1 & cumtoxC>=toxcon) break
if (newdose==1 & cumtoxC<toxcon){
if (cumtox>upperlimit) break
if (cumtox<=upperlimit)  doselevel<-newdose 
if (cumtox<lowerlimit)  doselevel<-newdose+1 
}
if (newdose==stopdose & cumtoxC>=toxcon) 
   {doselevel<-newdose-1
    stopdose<-newdose-1
}
if (newdose==stopdose & cumtoxC<toxcon){
if (cumtox>upperlimit & newdose>1) doselevel<-newdose-1
if (cumtox>=tootoxic)  stopdose<-newdose-1
if (cumtox<=upperlimit)   doselevel<-newdose 
}
newdose<-doselevel
if (cohortSize*length(currentlevel$tox)>=MTDss & cumtoxC<toxcon & cumtox>=lowerlimit) break
if (cohortSize*length(currentlevel$tox)>=MTDss & cumtox<lowerlimit & stopdose==doselevel) break
}
cumtox.m1<-cumtox.old[-c(1)]
cumtoxC.m1<-cumtoxC.old[-c(1)]
simNo<-rep(j, length(output$tox))
results<-cbind(simNo,output, cumtox.m1, cumtoxC.m1)
results.old<-rbind(results.old, results)
}
simresults<-data.frame(results.old[-c(1),],row.names=NULL)
simData<-(list(simresults=simresults, cohortSize=cohortSize,proba=proba,probc=probc, MTDss=MTDss,pTarget=pTarget,lowerlimit=lowerlimit,upperlimit=upperlimit,tootoxic=tootoxic, toxcon=toxcon, llactivity=llactivity, sim=sim))
DLdata<-teqr.DLdataconACT(simData=simData)
MTDdata<-teqr.MTDdataconACT(simData=simData, DLdata=DLdata)
OperCharact<-teqr.OperCharact(simData=simData, DLdata=DLdata,MTDdata=MTDdata)
OperCharact
}

#orders the clinical trial simulation results by doselevel
teqr.DLdata<-function(simData=simData){
sim<-simData$sim
SDR<-simData$simresults
re.old<-rep(NA,7)
for (i in 1:sim){
#print(i)
output<-SDR[SDR$simNo==i, ]
mindose<-min(output$doselevel)
maxdose<-max(output$doselevel)
for (j in mindose:maxdose){
outputd<-output[output$doselevel==j, ]
simNo<-outputd$simNo
doselevel<-outputd$doselevel
stox<-sum(outputd$tox)
dllength<-simData$cohortSize*length(outputd$tox)
if (dllength==0) next
toxl<-binom.test(stox,dllength)$conf[1]
toxu<-binom.test(stox,dllength)$conf[2]
toxest<-binom.test(stox,dllength)$estimate
re1<-cbind(simNo,doselevel,stox,dllength,toxl,toxu,toxest)
re<-re1[1,]
re.old<-rbind(re.old,re)
}
}
results<-data.frame(re.old[-c(1),], row.names=NULL)
names(results)<-c('simNo','doselevel','stox','dllength','toxl','toxu','toxest')
return(results=results)
}

#orders the clinical trial simulation results by doselevel
teqr.DLdataconACT<-function(simData=simData){
sim<-simData$sim
SDR<-simData$simresults
re.old<-rep(NA,9)
for (i in 1:sim){
#print(i)
output<-SDR[SDR$simNo==i, ]
mindose<-min(output$doselevel)
maxdose<-max(output$doselevel)
for (j in mindose:maxdose){
outputd<-output[output$doselevel==j, ]
simNo<-outputd$simNo
doselevel<-outputd$doselevel
stox<-sum(outputd$tox)
stoxC<-sum(outputd$toxC)
dllength<-simData$cohortSize*length(outputd$tox)
if (dllength==0) next
toxl<-binom.test(stox,dllength)$conf[1]
toxu<-binom.test(stox,dllength)$conf[2]
toxest<-binom.test(stox,dllength)$estimate
toxCest<-binom.test(stoxC,dllength)$estimate
re1<-cbind(simNo,doselevel,stox,stoxC,dllength,toxl,toxu,toxest,toxCest)
re<-re1[1,]
re.old<-rbind(re.old,re)
}
}
results<-data.frame(re.old[-c(1),], row.names=NULL)
names(results)<-c('simNo','doselevel','stox','stoxC','dllength','toxl','toxu','toxest', 'toxCest')
return(results=results)
}


#MTDdataset
teqr.MTDdatatox<-function(simData=simData,DLdata=DLdata){
doselevel=rep(NA,simData$sim)
for (i in 1:simData$sim){
output<-DLdata[DLdata$simNo==i, ]
iv<-rep(i,length(output$doselevel))
dat<-data.frame(cbind(c(iv),output$doselevel,round(c(pava(output$stox/output$dllength)),2)))
names(dat)<-c('i','dl','tox')
dat$tt<- ifelse(dat$tox<simData$tootoxic,1,0)
#print(dat)
if (sum(dat$tt)>0) {
newdat<-dat[dat$tt>0,] 
newdat$diff<-abs(newdat$tox-simData$pTarget)
newdat$mindiff<-min(newdat$diff)
mindiffdat<-newdat[newdat$diff==newdat$mindiff,]
doselevel[i]<-max(mindiffdat$dl)
}
if (sum(dat$tt)==0) doselevel[i]<-100
}
dose0<-DLdata$simNo[DLdata$doselevel==1 & DLdata$toxest>simData$upperlimit]
simNo<-seq(1:simData$sim)
maxpava<-data.frame(cbind(cbind(simNo),cbind(doselevel)))
for(i in dose0) maxpava$doselevel[i]<-0
#maxpava
Remax<-merge(maxpava,DLdata)
o<-order(Remax$simNo)
Remax.o<-Remax[o,]
return(Remax.o=Remax.o)
}


#MTDdataset
teqr.MTDdataconACT<-function(simData=simData,DLdata=DLdata){
doselevel=rep(NA,simData$sim)
for (i in 1:simData$sim){
output<-DLdata[DLdata$simNo==i, ]
output$tt<- ifelse(output$toxCest<simData$toxcon,1,0)
#print(output)
#removes simulations where all results are too toxic
if(sum(output$tt)<1)doselevel[i]<- -1
if(sum(output$tt)>0){
newoutput<-output[output$tt>0,]
iv<-rep(i,length(newoutput$doselevel))
newdat<-data.frame(cbind(c(iv),newoutput$doselevel,round(c(pava(newoutput$stox/newoutput$dllength)),2)),newoutput$toxCest)
names(newdat)<-c('i','dl','tox','toxC')
newdat$diff<-abs(newdat$tox-simData$pTarget)
newdat$mindiff<-min(newdat$diff)
mindiffdat<-newdat[newdat$diff==newdat$mindiff,]
#print(mindiffdat)
doselevel[i]<-min(mindiffdat$dl)
actdat<-mindiffdat[mindiffdat$dl==doselevel[i],]
#removes simulations were results are not active enough
if (actdat$tox<simData$llactivity) doselevel[i]<-100
}
}
simNo<-seq(1:simData$sim)
maxpava<-data.frame(cbind(cbind(simNo),cbind(doselevel)))
maxpava
Remax<-merge(maxpava,DLdata)
o<-order(Remax$simNo)
Remax.o<-Remax[o,]
return(Remax.o=Remax.o)
}

#Results for OperChartox
teqr.OperChartox<-function(simData=simData,DLdata=DLdata,MTDdata=MTDdata){
totalN<-rep(NA,simData$sim)
for (i in 1:simData$sim){ 
totalN[i]<-sum(DLdata$dllength[DLdata$simNo==i])
}
MedN<-median(totalN)
out<-table(simData$simresults$simNo, simData$simresults$doselevel)
NoPatients<-(simData$cohortSize*margin.table(out,2))/simData$sim
sstox<-rep(NA,simData$sim)
for (i in 1:simData$sim) sstox[i]<-sum(DLdata$stox[DLdata$simNo==i])
sdllength<-rep(NA,simData$sim)
for (i in 1:simData$sim) sdllength[i]<-sum(DLdata$dllength[DLdata$simNo==i])
DLTrate.dat<-data.frame(cbind(sstox,sdllength))
DLTrate<-sstox/sdllength
MeanDLTrate<-mean(DLTrate)
MeanToxRate<-mean(MTDdata$toxest)
MeanCIlength<-mean(MTDdata$toxu-MTDdata$toxl)
#percent of simulation achieving a defined samplesize
stopNo.m1=simData$MTDss-1
stopNo<-ifelse(MTDdata$dllength>stopNo.m1,1,0)
PropObd<-sum(stopNo)/simData$sim
NoMTD<-simData$sim-length(MTDdata$toxest)
Retab<-table(MTDdata$doselevel)
NoTrialsMTD<-prop.table(Retab)
NoTrialsMTD1<-prop.table(cbind(t(Retab),NoMTD))
names(simData$simresults)<-c("simNo","doselevel","tox","cumtox")
octox<-list(sim=simData$sim,MedN=MedN,NoPatients=NoPatients,MeanDLTrate=MeanDLTrate,NoTrialsMTD=NoTrialsMTD,NoTrialsMTD1=NoTrialsMTD1,MeanToxRate=MeanToxRate,MeanCIlength=MeanCIlength,PropObd=PropObd, NoMTD=NoMTD,simData=simData, DLdata=DLdata, MTDdata=MTDdata)
class(octox) <- "teqrOCtox"
    octox
}

#Results from OperCharact
teqr.OperCharact<-function(simData=simData,DLdata=DLdata,MTDdata=MTDdata){
totalN<-rep(NA,simData$sim)
for (i in 1:simData$sim){ 
totalN[i]<-sum(DLdata$dllength[DLdata$simNo==i])
}
MedN<-median(totalN)
out<-table(simData$simresults$simNo, simData$simresults$doselevel)
NoPatients<-(simData$cohortSize*margin.table(out,2))/simData$sim
sstox<-rep(NA,simData$sim)
for (i in 1:simData$sim) sstox[i]<-sum(DLdata$stox[DLdata$simNo==i])
sdllength<-rep(NA,simData$sim)
for (i in 1:simData$sim) sdllength[i]<-sum(DLdata$dllength[DLdata$simNo==i])
DLTrate.dat<-data.frame(cbind(sstox,sdllength))
DLTrate<-sstox/sdllength
MeanDLTrate<-mean(DLTrate)
MeanToxRate<-mean(MTDdata$toxest)
MeanCIlength<-mean(MTDdata$toxu-MTDdata$toxl)
MeanToxCestRate<-mean(MTDdata$toxCest)
#percent of simulation achieving a defined samplesize
stopNo.m1=simData$MTDss-1
stopNo<-ifelse(MTDdata$dllength>stopNo.m1,1,0)
PropObd<-sum(stopNo)/simData$sim
NoRP2D<-simData$sim-length(MTDdata$toxest)
Retab<-table(MTDdata$doselevel)
NoTrialsMTD<-prop.table(Retab)
NoTrialsMTD1<-prop.table(cbind(t(Retab),NoRP2D))
names(simData$simresults)<-c("simNo","doselevel","act","toxC","cumact","cumtoxC")
names(DLdata)<-c("simNo", "doselevel","sact","stoxC", "dllength", "actl","actu", "actest","toxCest")
names(MTDdata)<-c("simNo", "doselevel","sact","stoxC", "dllength", "actl","actu", "actest","toxCest")
ocact<-list(sim=simData$sim,MedN=MedN,NoPatients=NoPatients,MeanDLTrate=MeanDLTrate,NoTrialsMTD=NoTrialsMTD,NoTrialsMTD1=NoTrialsMTD1,MeanToxRate=MeanToxRate,MeanCIlength=MeanCIlength,MeanToxCestRate=MeanToxCestRate, PropObd=PropObd, NoRP2D=NoRP2D,simData=simData, DLdata=DLdata, RP2Ddata=MTDdata)
class(ocact) <- "teqrOCact"
    ocact
}




#dose escalation guidelines
teqrDG<-function(TotalN,pTarget,eq1,eq2,tootoxic){
Nplus1<-TotalN+1
upperlimit<-pTarget+eq2
lowerlimit<-pTarget-eq1
prob<-matrix(NA, ncol=TotalN, nrow=Nplus1, byrow=TRUE)
for (i in 1:Nplus1){
prob[i,]<-((i-1)/seq(1:TotalN))
}
row<-seq(1:Nplus1)-1
col<-seq(1:TotalN)
row.names(prob)<-row
colnames(prob)<-col
prob
prob2<-matrix(NA, ncol=TotalN, nrow=Nplus1, byrow=TRUE)
for (i in 1:Nplus1){
for (j in 1:TotalN) {
if (prob[i,j]<lowerlimit) prob2[i,j]<-'E'
if (prob[i,j]>= lowerlimit & prob[i,j] <= upperlimit) prob2[i,j]<-'S'
if (prob[i,j]>upperlimit) prob2[i,j]<-'D'
if (prob[i,j]>=tootoxic) prob2[i,j]<-'DU'
if (prob[i,j]>1.00)  prob2[i,j]<-'  '
}
} 
row<-seq(1:Nplus1)-1
col<-seq(1:TotalN)
row.names(prob2)<-row
colnames(prob2)<-col
prob2
dg<-(list(probTable=round(prob,3), DoseGuideTable=prob2))
class(dg) <- "teqrDG"
   dg
}


#print function
print.teqrOCtox<-function (x,...) 
{   cat("is an object of class teqrOCtox.   \n")
    cat("Ave. no of Patients studied at each dose level \n")
    print(round(x$NoPatients,2)) 
    cat("Rate dose level is chosed as the MTD \n")
    print(round(x$NoTrialsMTD1,2)) 
    cat("Median Study Sample Size:", x$MedN,"\n")
    cat("Mean Study DLT Rate:", round(x$MeanDLTrate,2),"\n")
    cat("Mean DLT Rate at the MTD:", round(x$MeanToxRate,2),"\n")
    cat("Mean 95% binomial confidence interval length at the MTD:", round(x$MeanCIlength,2),"\n") 
    cat("Proportion of trials with MTD dose level sample size","\n")
    cat("   at or above the desired number:", round(x$PropObd,2),"\n")
    cat("No of simulated trials that do not determine an MTD:", x$NoMTD,"\n")  
    cat("No of simulated trials:", x$sim,"\n")  
    cat("The following simulation data sets are also contained","\n")
    cat("   within the TEQR object.\n")
    cat("   The simulation level data: simData \n")
    cat("   The dose level data: DLdata \n")
    cat("   The MTD level data: MTDdata \n")
    invisible(NULL)
}

#print function
print.teqrOCact<-function (x,...) 
{   cat("is an object of class teqrOCact.   \n")  
    cat("Ave. no of Patients studied at each dose level \n")
    print(round(x$NoPatients,2))   
    cat("Rate dose level is chosed as the RP2D\n")
    print(round(x$NoTrialsMTD1,2))
    cat("Median Study Sample Size:", x$MedN,"\n")
    cat("Mean Study Event Rate:", round(x$MeanDLTrate,2),"\n")
    cat("Mean Event Rate at the RP2D:", round(x$MeanToxRate,2),"\n")
    cat("Mean 95% binomial confidence interval length at the RP2D:", round(x$MeanCIlength,2),"\n") 
    cat("Mean Toxicity Rate at the RP2D:", round(x$MeanToxCestRate,2),"\n")
    cat("Proportion of trials with RP2D dose level sample size","\n")
    cat("   at or above the desired number:", round(x$PropObd,2),"\n")
    cat("No of simulated trials that do not determine a RP2D:", x$NoRP2D,"\n")  
    cat("No of simulated trials:", x$sim,"\n")  
    cat("The following simulation data sets are also contained","\n")
    cat("   within the TEQR object.\n")
    cat("   The simulation level data: simData \n")
    cat("   The dose level data: DLdata \n")
    cat("   The RP2D level data: RP2Ddata \n")
    invisible(NULL)
}

#print function
print.teqrDG<-function (x,...) 
{   cat("is an object of class teqrDG which contains the dose escalation/de-escalation\n") 
    cat("guidelines table. The rows represent number of subjects that have\n")
    cat("experienced a DLT/event and the columns represent the number of subjects on\n")
    cat("the current dose level. The letter codes are the guidelines and the letters\n")
    cat("are defined as follows; E - escalate, S - Stay, D - De-escalate,\n")  
    cat("and DU - De-escalate and do not return to this dose. \n")  
    print(x$DoseGuideTable) 
    cat("                              \n")      
    cat("To print the dosing escalation/de-escalation guidelines\n")
    cat("or the underlying toxicity/event probabilities as objects\n") 
    cat("type x$DoseGuideTable and x$probTable, respectively.\n")
    invisible(NULL)
}



