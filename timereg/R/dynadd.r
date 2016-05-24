dynregBase<-function(times,status,response,fdata,designX,designA,id, clusters,
bhat=NULL,b=1,sim=1,antsim=1000,resample=0,smoothXX=0,
weighted.test=0)
{
Ntimes <- length(times)
designA<-as.matrix(designA); 
pa <- as.integer(dim(designA)[2]); na <- as.integer(dim(designA)[1])
designX<-as.matrix(designX); 
px <- as.integer(dim(designX)[2]); nx <- as.integer(dim(designX)[1])

#if (nx!=na) print(" A design og B designs er ikke ens\n");
w<-rep(1,nx); mw<-0;  
vcum.ly<-cum.ly<-cum0<-cumf<-cum.ms<-vcum0<-vcumf<-robvcumf<-matrix(0,Ntimes,px+1); 
ly<-NULL;

if (resample==1) cumBi<-matrix(0,Ntimes,fdata$antpers*px) else cumBi<-0;
test<-matrix(0,antsim,3*px); testOBS<-rep(0,3*px); testval<-c();
unifCI<-c();
rani<--round(runif(1)*10000)

# 50 tilfaeldige score processer til test H: b(t)=b returneres
if (sim==1) simUt<-matrix(0,Ntimes,50*px) else simUt<-NULL;
Ut<-matrix(0,Ntimes,px+1);

if (is.null(id) == TRUE) {
   antpers <- length(time); id <- 0:(antpers - 1);
} else { pers <- unique(id); antpers <- length(pers);
    id<-as.integer(factor(id, labels = 1:(antpers))) - 1; 
}

   clusters<-id; antclust<-antpers;
   fdata$antpers <- antpers; fdata$antclust <- antclust

out.aalen<-aalenBaseC(times,fdata,designA,status,id,clusters); 

###aalenBase(times,fdata,designX,status,id,clusters,robust=0,sim=0,retur=0,antsim=1000,weighted.test=1,
###	  covariance=0,resample.iid=0,namesX=NULL,silent=0,weights=NULL,entry=NULL,offsets=0)

if (!is.null(bhat)) xval<-bhat[,1] else
xval<-seq(times[2],times[Ntimes],length=30); 
smooth.aalen<-CsmoothB(out.aalen$cum,xval,b); 
nxval<-length(xval); 
#print(smooth.aalen[,1:4]); print(bhat); 

if (is.null(bhat)) { bhat<-matrix(0,nxval,px+1);  bhat[,1]<-smooth.aalen[,1];}
if (length(b)==1) b<-rep(b,nxval); 
bhatny<-bhat

nparout<-.C("dynadd",
as.double(times),as.double(response),as.integer(Ntimes),
as.double(designX),as.integer(nx),as.integer(px),
as.double(designA),as.integer(na),as.integer(pa),
as.double(smooth.aalen),as.double(bhat),as.double(bhatny),
as.integer(nxval),as.integer(fdata$antpers),as.double(fdata$start),
as.double(fdata$stop),
as.double(cum0), as.double(cumf), as.double(cum.ms),
as.double(vcum0), as.double(vcumf), as.double(robvcumf), 
as.double(w),as.integer(mw),as.integer(rani),
as.integer(sim),as.integer(antsim),
as.double(cumBi),
as.double(test),as.double(testOBS),
as.integer(status),as.double(Ut),as.double(simUt),
as.double(b),as.double(cum.ly),as.integer(resample),as.integer(id),
as.integer(smoothXX),as.integer(weighted.test),as.double(vcum.ly),
as.integer(clusters),as.integer(fdata$antclust)) ### , PACKAGE="timereg")

cum0 <-matrix(nparout[[17]],Ntimes,px+1);
cumf <-matrix(nparout[[18]],Ntimes,px+1);
cum.ms <-matrix(nparout[[19]],Ntimes,px+1);
vcum0 <-matrix(nparout[[20]],Ntimes,px+1);
vcumf <-matrix(nparout[[21]],Ntimes,px+1); 
robvcumf <-matrix(nparout[[22]],Ntimes,px+1);
cum.ly <-matrix(nparout[[35]],Ntimes,px+1) 
vcum.ly <-matrix(nparout[[40]],Ntimes,px+1) 

if (resample==1) {
cumBi<-matrix(nparout[[28]],Ntimes,fdata$antpers*px);
cumBI<-list();
for (i in (0:(fdata$antpers-1))*px)
cumBI[[i/px+1]]<-as.matrix(cumBi[,i+(1:px)]);  } else cumBI<-NULL;

if (sim==1) {
Uit<-matrix(nparout[[33]],Ntimes,50*px);
UIt<-list(); for (i in (0:49)*px) UIt[[i/px+1]]<-Uit[,i+(1:px)];
Ut<-matrix(nparout[[32]],Ntimes,(px+1));
test<-matrix(nparout[[29]],antsim,3*px);
testOBS<-nparout[[30]];
for (i in 1:(3*px)) testval<-c(testval,pval(test[,i],testOBS[i]))
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
} else {test<-NULL; unifCI<-NULL; Ut<-NULL; UIt<-NULL;
pval.testBeqC.is<-NULL; obs.testBeqC.is<-NULL; sim.testBeqC.is<-NULL;
pval.testBeq0<-NULL;pval.testBeqC<-NULL;
obs.testBeq0<-NULL;obs.testBeqC<-NULL;
sim.testBeq0<-NULL;sim.testBeqC<-NULL; }

out <- list(cum=cumf,var.cum=vcumf,robvar.cum=robvcumf,
cum0=cum0,var.cum0=vcum0,
cum.ms=cum.ms,var.cum.ms=vcumf,
cum.ly=cum.ly,var.cum.ly=vcum.ly,
B.iid=cumBI,
pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
pval.testBeqC.is=pval.testBeqC.is,
obs.testBeqC.is=obs.testBeqC.is,
sim.testBeqC.is=sim.testBeqC.is,
obs.testBeq0=obs.testBeq0,obs.testBeqC=obs.testBeqC,
sim.testBeq0= sim.testBeq0,sim.testBeqC=sim.testBeqC,
conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt)
return(out)
}

semiregBase<-function(times,status,response,fdata,designX,
designG,designA,id,clusters,gamma=0,bhat=NULL,b=1,sim=1,antsim=1000,
resample=0,weighted.test=0)
{
Ntimes <- length(times); maxtime<-times[Ntimes]; 
designX<-as.matrix(designX); designG<-as.matrix(designG); 
designA<-as.matrix(designA); 
px <- as.integer(dim(designX)[2]); nx <- as.integer(dim(designX)[1])
pg <- as.integer(dim(designG)[2]); ng <- as.integer(dim(designG)[1])
pa <- as.integer(dim(designA)[2]); nar <- as.integer(dim(designA)[1])

if (nx!=nar) print("Aalen-design and X designs are not consistent\n");
if (nx!=ng) print("Semi-design and and X designs are not consistent\n");

if (resample==1) {
gamma.iid<-matrix(0,fdata$antpers,pg); 
B.iid<-matrix(0,Ntimes,fdata$antpers*px) }
else { B.iid<-gamma.iid<-NULL; }

w<-rep(1,nx); mw<-0;  
cumly<-cum0<-cumf<-cumMS<-vcum0<-vcumf<-robvcumf<-matrix(0,Ntimes,px+1); 
test<-matrix(0,antsim,3*px); testOBS<-rep(0,3*px); testval<-c();
unifCI<-c(); rani<--round(runif(1)*10000)

# 50 tilfaeldige score processer til test H: b(t)=b returneres
if (sim==1) simUt<-matrix(0,Ntimes,50*px) else simUt<-NULL;
Ut<-matrix(0,Ntimes,px+1); resid<-matrix(0,Ntimes,fdata$antpers*px) 

if (is.null(id) == TRUE) {
   antpers <- length(time); id <- 0:(antpers - 1);} else 
   { pers <- unique(id); antpers <- length(pers);
    id<-as.integer(factor(id, labels = 1:(antpers))) - 1; }

   clusters<-id; antclust<-antpers;
   fdata$antpers <- antpers; fdata$antclust <- antclust


out.aalen<-aalenBaseC(times,fdata,designA,status,id,clusters); 
xval<-seq(times[2],times[Ntimes],length=30); 
smooth.aalen<-CsmoothB(out.aalen$cum,xval,b); 
#print(smooth.aalen); 
naval<-nrow(smooth.aalen); 

if (length(b)==1) b<-rep(b,naval); 
if (is.null(bhat)) {bhat<-matrix(0,naval,px+pg+1); 
                    bhat[,1]<-smooth.aalen[,1];}

nxval<-nrow(bhat); bhatny<-bhat; 
if (sum(gamma)==0) gamma<-apply(bhat[,(px+2):(px+pg+1)],1,mean);

#print("Prelim estimate of gamma based on non-parametric model is"); 
gamma2<-gamly<-gamkor<-gameffi<-gameffims<-gamma
cum0<-cumf<-cumef<-rvcum<-rvcumef<-matrix(0,Ntimes,px+1); 
Vgam0<-Vgamef<-Vkorgam<-robvargam<-robvargame<-VgammaLY<-Vgamma<-Vkorgam<-C<-matrix(0,pg,pg); 

semiout<-.C("semidynadd",
as.double(times),as.double(response),as.integer(Ntimes),
as.double(designX),as.integer(nx),as.integer(px),
as.double(designG),as.integer(ng),as.integer(pg),
as.double(designA),as.integer(nar),as.integer(pa),
as.double(smooth.aalen),as.integer(naval),as.double(bhat),
as.integer(nxval),as.integer(fdata$antpers),as.double(fdata$start),
as.double(fdata$stop),as.double(cum0),as.double(cumf),
as.double(cumMS), as.double(rvcum), as.double(rvcumef),
as.double(gamma),as.double(gamma2),as.double(gamly),
as.double(gamkor),as.double(gameffi),as.double(gameffims),
as.double(Vgamma),as.double(Vkorgam),as.double(Vgamef),
as.double(robvargam),as.double(robvargame), as.double(w),
as.integer(mw), as.integer(rani), as.integer(sim),
as.integer(antsim), as.double(resid),as.double(test),
as.double(testOBS), as.double(Ut),as.double(simUt),
as.double(b),as.integer(id), as.integer(status),
as.integer(weighted.test),as.double(VgammaLY),as.integer(clusters),
as.integer(fdata$antclust), as.integer(resample),as.double(gamma.iid),
as.double(B.iid),
PACKAGE="timereg") 

if (resample==1)  {
gamma.iid<-matrix(semiout[[54]],fdata$antclust,pg);
covit<-matrix(semiout[[55]],Ntimes,fdata$antclust*px); 
B.iid<-list(); 
for (i in (0:(fdata$antclust-1))*px) {
B.iid[[(i/px)+1]]<-as.matrix(covit[,i+(1:px)]);
###colnames(B.iid[[i/px+1]])<-namesX; }
###colnames(gamma.iid)<-namesZ
}
}

cum0 <-matrix(semiout[[20]],Ntimes,px+1); 
cumf <-matrix(semiout[[21]],Ntimes,px+1);
cum.ms <-matrix(semiout[[22]],Ntimes,px+1);
vcumf<-matrix(semiout[[23]],Ntimes,px+1);
robvcumf<-matrix(semiout[[23]],Ntimes,px+1);
vcum0 <-NULL; 
#cum.ly<-matrix(semiout[[36]],Ntimes,px+1);
cum.ly<-NULL; var.cum.ly<-NULL;

gamma<-matrix(semiout[[25]],pg,1); 
gamma2<-matrix(semiout[[26]],pg,1); 
gamly<-matrix(semiout[[27]],pg,1); 
gamkor<-matrix(semiout[[28]],pg,1); 
gameffi<-matrix(semiout[[29]],pg,1); 
gameffims<-matrix(semiout[[30]],pg,1);

Vgamma<-matrix(semiout[[31]],pg,pg); 
Vkorgam<-matrix(semiout[[32]],pg,pg); 
Vgamef<-matrix(semiout[[33]],pg,pg); 
robvargam<-matrix(semiout[[34]],pg,pg);
robvargame<-matrix(semiout[[35]],pg,pg);
VgammaLY<-matrix(semiout[[50]],pg,pg); 

if (resample==1) {
cumBi<-matrix(semiout[[41]],Ntimes,fdata$antpers*px); 
cumBI<-list();
for (i in (0:(fdata$antpers-1))*px)
cumBI[[i/px+1]]<-cumBi[,i+(1:px)] } else cumBI<-NULL;

if (sim==1) {
test<-matrix(semiout[[42]],antsim,3*px); testOBS<-semiout[[43]];
Ut<-matrix(semiout[[44]],Ntimes,(px+1));
Uit<-matrix(semiout[[45]],Ntimes,50*px);
UIt<-list(); for (i in (0:49)*px) UIt[[i/px+1]]<-Uit[,i+(1:px)];
for (i in 1:(3*px)) testval<-c(testval,pval(test[,i],testOBS[i]))
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
} else {test<-NULL; unifCI<-NULL; Ut<-NULL; UIt<-NULL;
pval.testBeq0<-NULL; obs.testBeq0<-NULL; sim.testBeq0<-NULL;
pval.testBeqC<-NULL; obs.testBeqC<-NULL; sim.testBeqC<-NULL; 
pval.testBeqC.is<-NULL; obs.testBeqC.is<-NULL; sim.testBeqC.is<-NULL;
}

#gamma.ef=gameffi,gamma.efms=gameffims,
#var.gamma.ef=Vgamef,robvar.gamma.ef=robvargame, 

out <- list(cum=cumf,var.cum=vcumf,robvar.cum=robvcumf,cum.ms=cum.ms,
cum0=cum0,var.cum0=vcum0,cum.ly=cum.ly,var.cum.ly=NULL,
gamma0=gamma,var.gamma0=Vgamma,gamma.ly=gamly,var.gamma.ly=VgammaLY,
gamma=gamkor,gamma.ms=gamma2,var.gamma=Vkorgam,var.gamma.ms=Vkorgam,
robvar.gamma=robvargam,
residuals=cumBI,
pval.testBeq0=pval.testBeq0, pval.testBeqC=pval.testBeqC,
pval.testBeqC.is=pval.testBeqC.is, obs.testBeqC.is=obs.testBeqC.is,
sim.testBeqC.is=sim.testBeqC.is,
obs.testBeq0=obs.testBeq0,obs.testBeqC=obs.testBeqC,
sim.testBeq0= sim.testBeq0,sim.testBeqC=sim.testBeqC,
conf.band=unifCI,test.procBeqC=Ut,sim.test.procBeqC=UIt,
gamma.iid=gamma.iid, B.iid=B.iid)
return(out)
}
