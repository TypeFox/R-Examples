timecoxBase<-function(times,fdata,designX,status,id,bhat,
method="basic",band=1,degree=1,it=20,sim=1,retur=0,robust=1,
antsim=1000,sim2=sim2,weighted.test=0,covariance=0)
{
Ntimes <- length(times); designX<-as.matrix(designX); 
if(is.matrix(designX) == TRUE) p<-as.integer(dim(designX)[2]);
if(is.matrix(designX) == TRUE) nx<-as.integer(dim(designX)[1]);
bhat<-as.matrix(bhat); schoen<-FALSE; 
nb<-as.integer(dim(bhat)[1]);
cum<-matrix(0,Ntimes,p+1); Vcum<-Vcum2<-rvcu<-cum; 
nullresid<-(-1);
if (method=="breslow") pdim<-p+1 else pdim<-p; 
band<-matrix(band,nb,p+1); 

if (retur==1) {cumAi<-matrix(0,Ntimes,fdata$antpers*1) 
               cumAiid<-matrix(0,Ntimes,fdata$antpers*1); } else { 
	       cumAi<-0; cumAiid<-0; }

if (covariance==1) covs<-matrix(0,Ntimes,p*p) else covs<-0;



if (sim>=1) {test<-matrix(0,5*pdim,antsim); testOBS<-rep(0,5*pdim);}
else {test<-0; testOBS<-0;} 
rani<--round(runif(1)*10000);  testval<-c();  unifCI<-c(); 

# 50 tilfaeldige score processer til test H: b(t)=b returneres
if (sim>=1) simUt<-matrix(0,Ntimes,50*pdim) else simUt<-NULL;
Ut<-matrix(0,Ntimes,pdim+1);

if (method=="basic") {
timeout<- .C("OStimecox",
as.double(times), as.integer(Ntimes),as.double(designX), as.integer(nx),as.integer(p),as.integer(fdata$antpers), 
as.double(fdata$start),as.double(fdata$stop),as.integer(nb), as.double(bhat), as.double(cum),as.double(Vcum),
as.integer(it),as.double(band),as.integer(degree), as.integer(id), as.integer(status), as.integer(sim),
as.integer(antsim), as.double(cumAi), as.double(test), as.integer(rani),as.double(testOBS),as.double(Ut),
as.double(simUt), as.double(rvcu), as.integer(retur), as.integer(weighted.test),as.double(cumAiid), as.integer(robust),
as.integer(covariance),as.double(covs), PACKAGE="timereg");

if (covariance==1)  {
covit<-matrix(timeout[[32]],Ntimes,p*p);
cov.list<-list();
for (i in 1:Ntimes) cov.list[[i]]<-matrix(covit[i,],p,p);
} else cov.list<-NULL;


schoen<- 
obs.testBeqC.is1<- obs.testBeqC.is2<- 
pval.testBeqC.is1<- pval.testBeqC.is2<- 
sim.testBeqC.is1<- sim.testBeqC.is2<-NULL;

rvcu<-matrix(timeout[[26]],Ntimes,p+1);

if (retur==1) {
cumAi<-matrix(timeout[[20]],Ntimes,fdata$antpers*1);
cumAiid<-matrix(timeout[[29]],Ntimes,fdata$antpers*1);
cumAi<-list(time=times,dM=cumAi,dM.iid=cumAiid);} else cumAi<-NULL;

if (sim>=1) {
Ut<-matrix(timeout[[24]],Ntimes,(p+1));
Uit<-matrix(timeout[[25]],Ntimes,50*p);
UIt<-list(); for (i in (0:49)*p) UIt[[i/p+1]]<-as.matrix(Uit[,i+(1:p)]);
test<-matrix(timeout[[21]],antsim,5*p); test<-test[,1:(3*p)]; 
testOBS<-timeout[[23]];
for (i in 1:(3*p)) testval<-c(testval,pval(test[,i],testOBS[i]))
for (i in 1:p)
unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
pval.testBeq0<-as.vector(testval[1:p]);
pval.testBeqC<-as.vector(testval[(p+1):(2*p)]);
pval.testBeqC.is<-as.vector(testval[(2*p+1):(3*p)]);
obs.testBeq0<-as.vector(testOBS[1:p]);
obs.testBeqC<-as.vector(testOBS[(p+1):(2*p)]);
obs.testBeqC.is<-as.vector(testOBS[(2*p+1):(3*p)]);
sim.testBeq0<-as.matrix(test[,1:p]);
sim.testBeqC<-as.matrix(test[,(p+1):(2*p)]);
sim.testBeqC.is<-as.matrix(test[,(2*p+1):(3*p)]);
} else {
test<- unifCI<- Ut<- UIt<-
pval.testBeq0<- 
obs.testBeq0<- sim.testBeq0<-
pval.testBeqC<- obs.testBeqC<- sim.testBeqC<- 
obs.testBeqC.is<- pval.testBeqC.is<- 
sim.testBeqC.is<- test.procBeqC<- test.procBeqC.is<-NULL; 
}
}
else if (method=="breslow") {
schoen<-matrix(0,Ntimes,p+1); 
cumlam<-cum<-Vcum<-rvcu<-matrix(0,Ntimes,p+2); 

timeout<-.C("OSbreslow", 
as.double(times),as.integer(Ntimes),as.double(designX),
as.integer(nx),as.integer(p),as.integer(fdata$antpers), 
as.double(fdata$start),as.double(fdata$stop),as.integer(nb),
as.double(bhat),as.double(cum),as.double(Vcum),
as.integer(it),as.double(band),as.integer(degree),
as.double(schoen),as.integer(sim),as.integer(antsim),
as.double(test),as.integer(rani),as.double(testOBS),
as.double(rvcu),as.double(cumlam),as.integer(nullresid),
as.integer(status),as.integer(id),as.integer(sim2),
as.double(Ut),as.double(simUt),as.integer(weighted.test),
as.integer(robust),
PACKAGE="timereg"); 

cov.list<-NULL; 
schoen<-matrix(timeout[[16]],Ntimes,p+1); 
rvcu<-matrix(timeout[[22]],Ntimes,p+2);

if (sim2!=1) {
obs.testBeqC.is1<- obs.testBeqC.is2<- 
pval.testBeqC.is1<- pval.testBeqC.is2<- 
sim.testBeqC.is1<- sim.testBeqC.is2<-NULL;
}

cumAi<-NULL;  # calculated for time basic version

if (sim>=1) {
Ut<-matrix(timeout[[28]],Ntimes,(pdim+1));
Uit<-matrix(timeout[[29]],Ntimes,50*pdim);
UIt<-list(); for (i in (0:49)*pdim) UIt[[i/pdim+1]]<-as.matrix(Uit[,i+(1:pdim)]);
test<-matrix(timeout[[19]],antsim,5*pdim);
testOBS<-timeout[[21]]; testval<-c(); unifCI<-c(); 
for (i in 1:(5*pdim)) testval<-c(testval,pval(test[,i],testOBS[i]))
for (i in 1:pdim) unifCI<-c(unifCI,percen(test[,i],0.95)); 
pval.testBeq0<-as.vector(testval[1:pdim]);
pval.testBeqC<-as.vector(testval[(pdim+1):(2*pdim)]);
pval.testBeqC.is<-as.vector(testval[(2*pdim+1):(3*pdim)]);
obs.testBeq0<-as.vector(testOBS[1:pdim]);
obs.testBeqC<-as.vector(testOBS[(pdim+1):(2*pdim)]);
obs.testBeqC.is<-as.vector(testOBS[(2*pdim+1):(3*pdim)]);
sim.testBeq0<-as.matrix(test[,1:pdim]);
sim.testBeqC<-as.matrix(test[,(pdim+1):(2*pdim)]);
sim.testBeqC.is<-as.matrix(test[,(2*pdim+1):(3*pdim)]);
sim2<-0
if (sim2==1) {
pval.testBeqC.is1<-as.vector(testval[(3*pdim+1):(4*pdim)]);
pval.testBeqC.is2<-as.vector(testval[(4*pdim+1):(5*pdim)]);
obs.testBeqC.is1<-as.vector(testOBS[(3*pdim+1):(4*pdim)]);
obs.testBeqC.is2<-as.vector(testOBS[(4*pdim+1):(5*pdim)]);
sim.testBeqC.is1<-as.matrix(test[,(3*pdim+1):(4*pdim)]);
sim.testBeqC.is2<-as.matrix(test[,(3*pdim+1):(4*pdim)]); }
} else {test<- unifCI<- Ut<- UIt<-
pval.testBeq0<-pval.testBeqC<-
obs.testBeq0<-obs.testBeqC<-
sim.testBeq0<-sim.testBeqC<- 
sim.testBeqC.is<- sim.testBeqC.is1<- sim.testBeqC.is2<- 
pval.testBeqC.is<- pval.testBeqC.is1<- pval.testBeqC.is2<- 
obs.testBeqC.is<- obs.testBeqC.is1<- obs.testBeqC.is2<-NULL; 
}
} else  stop("Methods are : breslow and basic\n"); 

bhat<-matrix(timeout[[10]],nb,pdim+1); 
cum<-matrix(timeout[[11]],Ntimes,pdim+1); 
Vcum<-matrix(timeout[[12]],Ntimes,pdim+1);

# additional output ? band=band,method=method,degree=degree,it=it,
# beta.t=bhat,
ud<-list(cum=cum,var.cum=Vcum,robvar.cum=rvcu,
residuals=cumAi,schoenfeld.residuals=schoen,
schoenfeld.residual=schoen,
obs.testBeq0=obs.testBeq0, pval.testBeq0=pval.testBeq0,
pval.testBeqC=pval.testBeqC, obs.testBeqC=obs.testBeqC,
obs.testBeqC.is=obs.testBeqC.is, pval.testBeqC.is=pval.testBeqC.is,
obs.testBeqC.is1=obs.testBeqC.is1,pval.testBeqC.is1=pval.testBeqC.is1,
obs.testBeqC.is2=obs.testBeqC.is2,pval.testBeqC.is2=pval.testBeqC.is2,
sim.testBeq0= sim.testBeq0,sim.testBeqC=sim.testBeqC,
conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt,
sim.testBeqC.is=sim.testBeqC.is, sim.testBeqC.is=sim.testBeqC.is1,
sim.testBeqC.is=sim.testBeqC.is2,covariance=cov.list)
return(ud); 
}

semicox<-function(times,fdata,designX,designG,status,id,bhat,gamma=0,
band=1,degree=1,it=10,method="basic",sim=0,retur=0,antsim=1000,
robust=1,weighted.test=0,covariance=0) 
{
Nalltimes <- length(times); 
Ntimes<-sum(status[(fdata$stop>times[1]) & (fdata$stop<=times[Nalltimes])])+1;
designX<-as.matrix(designX); designG<-as.matrix(designG); 
if(is.matrix(designX) == TRUE) px <- as.integer(dim(designX)[2])
if(is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])
if(is.matrix(designG) == TRUE) pg <- as.integer(dim(designG)[2])
if(is.matrix(designG) == TRUE) ng <- as.integer(dim(designG)[1])
bhat<-as.matrix(bhat); 
nb<-as.integer(dim(bhat)[1]);
band<-matrix(band,nb,as.integer(dim(bhat)[2])-1); 
if (method=="breslow") pdim<-px+1 else pdim<-px; 

if (covariance==1) covs<-matrix(0,Ntimes,px*px) else
covs<-0;


if (nx!=ng) print(" A design og B designs er ikke ens\n");
cum<-matrix(0,Ntimes,px+1); rvcu<-Vcum<-cum; 
if (sum(abs(gamma))==0) gamma<-rep(0,pg)  else gamma<-gamma; 
RobVargam<-Vargam<-matrix(0,pg,pg); 

test<-matrix(0,antsim,3*pdim); testOBS<-rep(0,3*pdim);
testval<-c(); rani<--round(runif(1)*10000);

# 50 tilfaeldige score processer til test H: b(t)=b returneres
if (sim>=1) simUt<-matrix(0,Ntimes,50*pdim) else simUt<-NULL;
Ut<-matrix(0,Ntimes,pdim+1);
if (retur==1) cumAi<-matrix(0,Ntimes,fdata$antpers) else cumAi<-0;
sim.testBeqC.is1<- sim.testBeqC.is2<-
obs.testBeqC.is1<- obs.testBeqC.is2<-
pval.testBeqC.is1<- pval.testBeqC.is2<-NULL;

# LOAD NECESSARY ROUTINES
#if (method=="basic" && system=="unix") dyn.load("timecox.so")
#if (method=="basic" && system=="linux") dyn.load("lintimecox.so")
#else if (method=="breslow" && system=="unix") dyn.load("breslow.so") 
#else if (method=="breslow" && system=="linux") dyn.load("linbreslow.so")
#else if (method=="breslow" && system=="windows") dyn.load("breslow.dll"); 


if (method=="basic") 
{
semiout<-.C("OSsemicox",
as.double(times),as.integer(Ntimes),
as.double(designX),
as.integer(nx),as.integer(px),as.double(designG),
as.integer(ng),as.integer(pg),as.integer(fdata$antpers),
as.double(fdata$start),as.double(fdata$stop),as.integer(nb),
as.double(bhat),as.double(cum),as.double(Vcum),
as.double(gamma),as.double(Vargam),as.double(band),
as.integer(degree), as.integer(it), as.double(RobVargam), 
as.double(rvcu), as.integer(sim), as.integer(antsim), 
as.integer(retur), as.double(cumAi), as.double(test), 
as.integer(rani), as.double(testOBS), as.integer(status), 
as.double(Ut), as.double(simUt), as.integer(id),
as.integer(weighted.test),as.integer(robust),
as.integer(covariance),as.double(covs),PACKAGE="timereg");

if (covariance==1)  {
covit<-matrix(semiout[[37]],Ntimes,px*px);
cov.list<-list();
for (i in 1:Ntimes) cov.list[[i]]<-matrix(covit[i,],px,px);
} else cov.list<-NULL;

bhat<-matrix(semiout[[13]],nb,px+1);
cum<-matrix(semiout[[14]],Ntimes,px+1);
Vcum<-matrix(semiout[[15]],Ntimes,px+1);
rvcu<-matrix(semiout[[22]],Ntimes,px+1);

gamma<-matrix(semiout[[16]],pg,1);
Vargam<-matrix(semiout[[17]],pg,pg);
RobVargam<-matrix(semiout[[21]],pg,pg);

if (retur==1) {
cumAi<-matrix(semiout[[26]],Ntimes,fdata$antpers);
cumAi<-list(time=times,dM=cumAi); 
#cumAI<-list();
#for (i in (0:(fdata$antpers-1))*p) cumAI[[i/p+1]]<-cumAi[,i+(1:p)]
} else cumAi<-NULL;


unifCI<-c();
if (sim>=1) {
test<-matrix(semiout[[27]],antsim,3*px);
testOBS<-semiout[[29]];
for (i in 1:(3*px))
testval<-c(testval,pval(test[,i],testOBS[i]))
Ut<-matrix(semiout[[31]],Ntimes,(px+1));
Uit<-matrix(semiout[[32]],Ntimes,50*px);
UIt<-list(); for (i in (0:49)*px) UIt[[i/px+1]]<-as.matrix(Uit[,i+(1:px)]);
for (i in 1:px)
unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
obs.testBeq0<-as.vector(testOBS[1:px]);
pval.testBeq0<-as.vector(testval[1:px]);
sim.testBeq0<-as.matrix(test[,1:px]);
obs.testBeqC<-as.vector(testOBS[(px+1):(2*px)]);
pval.testBeqC<-as.vector(testval[(px+1):(2*px)]);
sim.testBeqC<-as.matrix(test[,(px+1):(2*px)]);
pval.testBeqC.is<-as.vector(testval[(2*px+1):(3*px)]);
obs.testBeqC.is<-as.vector(testOBS[(2*px+1):(3*px)]);
sim.testBeqC.is<-as.matrix(test[,(2*px+1):(3*px)]);
} else {
test<- unifCI<-UIt<- Ut<-
pval.testBeq0<-pval.testBeqC<-
obs.testBeq0<- obs.testBeqC<-
obs.testBeqC.is<- pval.testBeqC.is<-
sim.testBeq0<- sim.testBeqC<- sim.testBeqC.is<-NULL; 
}
}
else if (method=="breslow") 
{

schoen<-0; cum<-Vcum<-rvcu<-matrix(0,Ntimes,px+2); 
cumAi<-NULL; 

semiout<-.C("semibreslow",
as.double(times),as.integer(Ntimes),
as.double(designX),
as.integer(nx),as.integer(px),as.double(designG),
as.integer(ng),as.integer(pg),as.integer(fdata$antpers),
as.double(fdata$start),as.double(fdata$stop),as.integer(nb),
as.double(bhat),as.double(cum), as.double(Vcum), 
as.double(rvcu),as.double(gamma),as.double(Vargam),
as.double(RobVargam),as.double(band),as.integer(degree),
as.integer(it),as.integer(sim),as.integer(antsim),
as.double(test),as.integer(rani),as.double(testOBS),
as.integer(status),as.integer(id),as.double(schoen),
as.double(simUt),as.double(Ut),PACKAGE="timereg")  

bhat<-matrix(semiout[[13]],nb,px+2);
cum <-matrix(semiout[[14]],Ntimes,px+2);
Vcum <-matrix(semiout[[15]],Ntimes,px+2);
rvcu<-matrix(semiout[[16]],Ntimes,px+2);
gamma<-matrix(semiout[[17]],pg,1);
Vargam<-matrix(semiout[[18]],pg,pg);
RobVargam<-matrix(semiout[[19]],pg,pg);

schoen<-NULL; 
sim.testBeqC.is1<- sim.testBeqC.is2<-
obs.testBeqC.is1<- obs.testBeqC.is2<-
pval.testBeqC.is1<- pval.testBeqC.is2<-NULL;

if (sim>=1) {
unifCI<-c();
test<-matrix(semiout[[25]],antsim,3*pdim);
testOBS<-semiout[[27]];
for (i in 1:(3*pdim)) testval<-c(testval,pval(test[,i],testOBS[i]))
Uit<-matrix(semiout[[31]],Ntimes,50*pdim);
UIt<-list(); for (i in (0:49)*pdim)
UIt[[i/pdim+1]]<-as.matrix(Uit[,i+(1:pdim)]);
Ut<-matrix(semiout[[32]],Ntimes,(pdim+1));
for (i in 1:pdim) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
obs.testBeq0<-as.vector(testOBS[1:pdim]);
pval.testBeq0<-as.vector(testval[1:pdim]);
sim.testBeq0<-as.matrix(test[,1:pdim]);
pval.testBeqC<-as.vector(testval[(pdim+1):(2*pdim)]);
obs.testBeqC<-as.vector(testOBS[(pdim+1):(2*pdim)]);
sim.testBeqC<-as.matrix(test[,(pdim+1):(2*pdim)]);
pval.testBeqC.is<-as.vector(testval[(2*pdim+1):(3*pdim)]);
obs.testBeqC.is<-as.vector(testOBS[(2*pdim+1):(3*pdim)]);
sim.testBeqC.is<-as.matrix(test[,(2*pdim+1):(3*pdim)]);
} else {
test<-unifCI<-UIt<-Ut<-pval.testBeq0<-pval.testBeqC<-
obs.testBeq0<-obs.testBeqC<-obs.testBeqC.is<-Fpval.testBeqC.is<-
sim.testBeq0<-sim.testBeqC<-sim.testBeqC.is<-NULL; 
}
cov.list<-NULL; 

} else stop("Methods are : breslow and basic\n"); 

# method=method
ud<-list(cum=cum,var.cum=Vcum,robvar.cum=rvcu,
gamma=gamma,var.gamma=Vargam,robvar.gamma=RobVargam,
residuals=cumAi,
pval.testBeq0=pval.testBeq0, obs.testBeq0=obs.testBeq0,
pval.testBeqC=pval.testBeqC, obs.testBeqC=obs.testBeqC,
pval.testBeqC.is=pval.testBeqC.is, obs.testBeqC.is=obs.testBeqC.is,
sim.testBeq0= sim.testBeq0,sim.testBeqC=sim.testBeqC,
sim.testBeqC.is=sim.testBeqC.is,
conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt,
covariance=cov.list)
return(ud);
}

pval<-function(simt,Otest)
{
simt<-sort(simt);
p<-sum(Otest<simt)/length(simt);
return(p)
}

percen<-function(x,per)
{ n<-length(x); tag<-round(n*per)+1; ud<-sort(x)[tag];
return(ud) }

