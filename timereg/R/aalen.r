aalenBase<-function(times,fdata,designX,status,id,clusters,robust=0,sim=0,retur=0,antsim=1000,weighted.test=1,covariance=0,resample.iid=0,namesX=NULL,silent=0,weights=NULL,entry=NULL,offsets=0) 
{
## {{{ setup variables
Ntimes<-length(times); designX<-as.matrix(designX); 
if(is.matrix(designX) == TRUE) p <- as.integer(dim(designX)[2])
if(is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])

if (is.diag(  t(designX) %*% designX  )==TRUE) stratum <- 1 else stratum <- 0

if (robust==0 & sim>=1)  robust<-1;  
cumint<-matrix(0,Ntimes,p+1); Vcumint<-cumint; robVar<-Vcumint; 
cumAi<-0;
if (retur==1) cumAi<-matrix(0,Ntimes,fdata$antpers) 
if (retur==2) cumAi<-rep(0,fdata$antpers); 
test<-matrix(0,antsim,3*p); testOBS<-rep(0,3*p); testval<-c(); 
unifCI<-c(); 

# 49 random score processes for testing H: b(t)=b returned 
if (sim>=1) simUt<-matrix(0,Ntimes,50*p) else simUt<-NULL;
Ut<-matrix(0,Ntimes,p+1); 
if (covariance==1) covs<-matrix(0,Ntimes,p*p) else covs<-0;
if (resample.iid==1) {
B.iid<-matrix(0,Ntimes,fdata$antclust*p) } else B.iid<-NULL;
if (sum(offsets)==0) mof <- 0 else mof <- 1; 
## }}}

aalenout<- .C("robaalen", ## {{{
as.double(times), as.integer(Ntimes),as.double(designX),
as.integer(nx),as.integer(p),as.integer(fdata$antpers),
as.double(fdata$start),as.double(fdata$stop),as.double(cumint), # 3
as.double(Vcumint),as.double(robVar),as.integer(sim),
as.integer(antsim),as.integer(retur),as.double(cumAi),
as.double(test),as.double(testOBS),as.integer(status),          # 6
as.double(Ut),as.double(simUt),as.integer(id),
as.integer(weighted.test),as.integer(robust),as.integer(covariance),
as.double(covs),as.integer(resample.iid),as.double(B.iid),      # 9
as.integer(clusters),as.integer(fdata$antclust),
as.integer(silent),as.double(weights),as.integer(entry),      # 11
as.integer(mof),as.double(offsets),as.integer(stratum)
,PACKAGE="timereg")
## }}}

## {{{ handling output
cumint <-matrix(aalenout[[9]],Ntimes,p+1);
Vcumint<-matrix(aalenout[[10]],Ntimes,p+1);
robVar<-matrix(aalenout[[11]],Ntimes,p+1);

if (covariance==1)  {
covit<-matrix(aalenout[[25]],Ntimes,p*p);
cov.list<-list();
for (i in 1:Ntimes) cov.list[[i]]<-matrix(covit[i,],p,p);
} else cov.list<-NULL;

if (resample.iid==1)  {
covit<-matrix(aalenout[[27]],Ntimes,fdata$antclust*p); 
B.iid<-list(); 
for (i in (0:(fdata$antclust-1))*p) {
B.iid[[i/p+1]]<-as.matrix(covit[,i+(1:p)]); 
colnames(B.iid[[i/p+1]])<-namesX; } 
}

cumAi<-NULL; 
if (retur==1) {
cumAi<-matrix(aalenout[[15]],Ntimes,fdata$antpers*1); 
cumAi<-list(time=times,dM=cumAi,dM.iid=cumAi); 
}
if (retur==2)  {cumAi<-aalenout[[15]]; cumAi<-list(dM=cumAi); }

if (sim>=1) {
Uit<-matrix(aalenout[[20]],Ntimes,50*p); 
UIt<-list(); for (i in (0:49)*p) 
UIt[[i/p+1]]<-as.matrix(Uit[,i+(1:p)]);
Ut<-matrix(aalenout[[19]],Ntimes,(p+1)); 
test<-matrix(aalenout[[16]],antsim,3*p);
testOBS<-aalenout[[17]];
for (i in 1:(3*p)) testval<-c(testval,pval(test[,i],testOBS[i]))
for (i in 1:p) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
pval.testBeq0<-as.vector(testval[1:p]);
pval.testBeqC<-as.vector(testval[(p+1):(2*p)]);
pval.testBeqC.is<-as.vector(testval[(2*p+1):(3*p)]);
obs.testBeq0<-as.vector(testOBS[1:p]);
obs.testBeqC<-as.vector(testOBS[(p+1):(2*p)]);
obs.testBeqC.is<-as.vector(testOBS[(2*p+1):(3*p)]);
sim.testBeq0<-as.matrix(test[,1:p]);
sim.testBeqC<-as.matrix(test[,(p+1):(2*p)]);
sim.testBeqC.is<-as.matrix(test[,(2*p+1):(3*p)]);
} else {test<-NULL; unifCI<-NULL; Ut<-NULL; UIt<-NULL; 
pval.testBeq0<-NULL;pval.testBeqC<-NULL; obs.testBeq0<-NULL;obs.testBeqC<-NULL;
sim.testBeq0<-NULL;sim.testBeqC<-NULL; 
sim.testBeqC.is<-NULL; pval.testBeqC.is<-NULL; 
obs.testBeqC.is<-NULL; 
}
## }}}

list(cum=cumint,var.cum=Vcumint,robvar.cum=robVar,residuals=cumAi,
pval.testBeq0=pval.testBeq0, obs.testBeq0=obs.testBeq0,
pval.testBeqC=pval.testBeqC, pval.testBeqC.is=pval.testBeqC.is,
obs.testBeqC=obs.testBeqC,obs.testBeqC.is=obs.testBeqC.is,
sim.testBeq0= sim.testBeq0,
sim.testBeqC=sim.testBeqC,sim.testBeqC.is=sim.testBeqC.is,
conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt,
covariance=cov.list,B.iid=B.iid,stratum=stratum)
}

const <- function(x) x 

pval<-function(simt,Otest)
{
simt<-sort(simt);
p<-sum(Otest<simt)/length(simt);
return(p)
}

percen<-function(x,per)
{ n<-length(x); tag<-round(n*per)+1; ud<-sort(x)[tag]; return(ud) }

semiaalen<-function(times,fdata,designX,designG,status,id,clusters,bhat=0,gamma=NULL,robust=0,sim=0,antsim=1000,weighted.test=1,retur=0,covariance=0,
resample.iid=0,namesX=NULL,namesZ=NULL,deltaweight=1,
silent=0,weights=NULL,entry=NULL,offsets=0)
{
## {{{ setting up variables 
Nalltimes <- length(times);  
Ntimes<-sum(status[(fdata$stop>times[1]) & (fdata$stop<=times[Nalltimes])])+1;
#print(Ntimes); print(Nalltimes); 
#print(times);  print(status);  print(cbind(fdata$stop,fdata$start,status))
designX<-as.matrix(designX); designG<-as.matrix(designG); 
if(is.matrix(designX) == TRUE) px <- as.integer(dim(designX)[2])
if(is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])
if(is.matrix(designG) == TRUE) pg <- as.integer(dim(designG)[2])
if(is.matrix(designG) == TRUE) ng <- as.integer(dim(designG)[1])
nb<-1; 
if(is.matrix(bhat)==TRUE) nb<-as.integer(dim(bhat)[1]); 
if(is.matrix(bhat)==FALSE) bhat<-matrix(0,nb,px+1); 

if (is.diag(  t(designX) %*% designX  )==TRUE) stratum <- 1 else stratum <- 0
if (covariance==1) covs<-matrix(0,Ntimes,px*px) else covs<-0;

if (resample.iid==1) {
gamma.iid<-matrix(0,fdata$antclust,pg); 
B.iid<-matrix(0,Ntimes,fdata$antclust*px) }
else { B.iid<-gamma.iid<-NULL; }

if (retur==1) cumAi<-matrix(0,Ntimes,fdata$antpers) else cumAi<-0;
if (nx!=ng) print(" A design og B designs er ikke ens\n");
cum<-matrix(0,Ntimes,px+1); 
Vcum<-cum; robVcum<-cum; 
Vargam2 <- Vargam<-matrix(0,pg,pg); RobVargam<-matrix(0,pg,pg); 
intZHZ<-matrix(0,pg,pg); gamma2 <- intZHdN<-rep(0,pg); 

test<-matrix(0,antsim,3*px); testOBS<-rep(0,3*px); testval<-c(); 

# 50 tilfaeldige score processer til test H: b(t)=b returneres
if (sim==1) simUt<-matrix(0,Ntimes,50*px) else simUt<-NULL;
Ut<-matrix(0,Ntimes,px+1); 

if (is.null(gamma)) fix.gamma <- 0 else fix.gamma <- 1; 
if (fix.gamma==1) { 
	if (length(gamma)!=pg) gamma <- rep(gamma[1],pg); 
} else gamma <- rep(0,pg);  
if (sum(abs(gamma))==0) gamma<-rep(0,pg)  else gamma<-gamma; 

if (sum(offsets)==0) mof<-0 else mof<-1;
## }}}

semiout<-.C("semiaalen", ## {{{
as.double(times),as.integer(Nalltimes),as.integer(Ntimes), # 1
as.double(designX),as.integer(nx),as.integer(px), # 2
as.double(designG),as.integer(ng),as.integer(pg), # 3
as.integer(fdata$antpers),as.double(fdata$start),as.double(fdata$stop), # 4
as.integer(nb),as.double(bhat), as.double(cum),   # 5
as.double(Vcum),as.double(robVcum),as.double(gamma),   # 6
as.double(Vargam),as.double(RobVargam), as.integer(sim),   # 7
as.integer(antsim), as.double(test),as.double(testOBS),   # 8
as.integer(robust), as.integer(status),as.double(Ut),   # 9
as.double(simUt),as.integer(id), as.integer(weighted.test),   # 10 
as.double(cumAi),as.integer(retur),as.integer(covariance),   # 11
as.double(covs), as.integer(resample.iid),as.double(gamma.iid),   # 12
as.double(B.iid), as.integer(clusters),as.integer(fdata$antclust),   # 13
as.double(intZHZ),as.double(intZHdN),as.integer(deltaweight),        # 14
as.integer(silent),as.double(weights),as.integer(entry),
as.integer(fix.gamma),as.integer(mof),as.double(offsets),
as.double(gamma2),as.double(Vargam2)
,PACKAGE="timereg"); 
## }}}

## {{{ handling output
if (resample.iid==1)  {
gamma.iid<-matrix(semiout[[36]],fdata$antclust,pg);
covit<-matrix(semiout[[37]],Ntimes,fdata$antclust*px); 
B.iid<-list(); 
for (i in (0:(fdata$antclust-1))*px) {
B.iid[[(i/px)+1]]<-as.matrix(covit[,i+(1:px)]);
colnames(B.iid[[i/px+1]])<-namesX; }
colnames(gamma.iid)<-namesZ
}

if (covariance==1)  {
covit<-matrix(semiout[[34]],Ntimes,px*px);
cov.list<-list();
for (i in 1:Ntimes) cov.list[[i]]<-matrix(covit[i,],px,px);
} else cov.list<-NULL;

bhat<-matrix(semiout[[14]],nb,px+1); cum <-matrix(semiout[[15]],Ntimes,px+1); 
Vcum <-matrix(semiout[[16]],Ntimes,px+1); 
robVcum <-matrix(semiout[[17]],Ntimes,px+1); robVcum[,1] <-cum[,1]; 
if (fix.gamma==0) gamma<-matrix(semiout[[18]],pg,1); 
Vargam<-matrix(semiout[[19]],pg,pg); 
robVargam<-matrix(semiout[[20]],pg,pg); 
intZHZ<-matrix(semiout[[40]],pg,pg); intZHdN<-matrix(semiout[[41]],pg,1); 
Vargam2 <-matrix(semiout[[50]],pg,pg); gamma2<-matrix(semiout[[49]],pg,1); 

if (retur==1) {
cumAi<-matrix(semiout[[31]],Ntimes,fdata$antpers*1); 
cumAi<-list(time=Vcum[,1],dM=cumAi); 
#cumAI<-list(); 
#for (i in (0:(fdata$antclust-1))*p) cumAI[[i/p+1]]<-cumAi[,i+(1:p)]
} else cumAi<-NULL;

unifCI<-c();
if (sim>=1) {
test<-matrix(semiout[[23]],antsim,3*px);
testOBS<-semiout[[24]];
for (i in 1:(3*px)) testval<-c(testval,pval(test[,i],testOBS[i]))
Uit<-matrix(semiout[[28]],Ntimes,50*px); 
UIt<-list(); for (i in (0:49)*px) UIt[[i/px+1]]<-Uit[,i+(1:px)];
Ut<-matrix(semiout[[27]],Ntimes,(px+1)); 
for (i in 1:px) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
pval.testBeq0<-as.vector(testval[1:px]);
pval.testBeqC<-as.vector(testval[(px+1):(2*px)]);
pval.testBeqC.is<-as.vector(testval[(2*px+1):(3*px)]);
obs.testBeq0<-as.vector(testOBS[1:px]);
obs.testBeqC<-as.vector(testOBS[(px+1):(2*px)]);
obs.testBeqC.is<-as.vector(testOBS[(2*px+1):(3*px)]);
sim.testBeq0<-as.matrix(test[,1:px]);
sim.testBeqC<-as.matrix(test[,(px+1):(2*px)]);
sim.testBeqC.is<-as.matrix(test[,(2*px+1):(3*px)]);
} else {test<-NULL; unifCI<-NULL;UIt<-NULL; Ut<-NULL;
pval.testBeq0<-NULL;pval.testBeqC<-NULL; 
obs.testBeq0<-NULL;obs.testBeqC<-NULL;
sim.testBeq0<-NULL; sim.testBeqC<-NULL;
sim.testBeqC.is<-NULL; pval.testBeqC.is<-NULL; 
obs.testBeqC.is<-NULL; 
}

tau=max(cum[,1])
## }}}

ud<-list(cum=cum,var.cum=Vcum,robvar.cum=robVcum,
gamma=gamma,var.gamma=Vargam,robvar.gamma=robVargam,residuals=cumAi,
pval.testBeq0=pval.testBeq0, obs.testBeq0=obs.testBeq0, 
pval.testBeqC=pval.testBeqC,pval.testBeqC.is=pval.testBeqC.is,
obs.testBeqC=obs.testBeqC, obs.testBeqC.is=obs.testBeqC.is,
sim.testBeq0= sim.testBeq0,
sim.testBeqC=sim.testBeqC,sim.testBeqC.is=sim.testBeqC.is,
conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt,
covariance=cov.list,B.iid=B.iid,gamma.iid=gamma.iid,
intZHZ=intZHZ,intZHdN=intZHdN,stratum=stratum); ### ,gamma2=gamma2/tau,var.gamma2=Vargam2/tau^2)
return(ud);
}
