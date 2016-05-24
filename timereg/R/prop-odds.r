prop.odds<-function(formula,data=sys.parent(),beta=NULL,
Nit=20,detail=0,start.time=0,max.time=NULL,id=NULL,n.sim=500,weighted.test=0,
profile=1,sym=0,baselinevar=1,clusters=NULL,max.clust=1000,weights=NULL)
{ ## {{{ 

out <- prop.odds.subdist(formula,data=data,beta=beta,cause=1,
Nit=Nit,detail=detail,start.time=start.time,max.time=max.time,id=id,n.sim=n.sim,
weighted.test=weighted.test,
profile=profile,sym=sym,cens.model="po",clusters=clusters,max.clust=max.clust,baselinevar=1,weights=weights)

return(out); 
} ## }}} 

prop.odds.gam<-function(formula,data=sys.parent(),beta=NULL,
Nit=10,detail=0,start.time=0,max.time=NULL,id=NULL,n.sim=500,weighted.test=0,
profile=1,sym=0,baselinevar=1,clusters=NULL,max.clust=1000)
{ ## {{{ 
id.call<-id; call<-match.call(); residuals<-0;  
robust<-0; ratesim<-0; 
resample.iid <- 1 # profile<-0; 
m<-match.call(expand.dots = FALSE);
m$sym<-m$profile<-m$max.time<-m$start.time<-m$weighted.test<-m$n.sim<-
m$id<-m$Nit<-m$detail<-m$beta <- m$baselinevar<-m$clusters <- m$max.clust  <-  NULL
if (n.sim==0) sim<-0 else sim<-1; 
antsim<-n.sim; 

Terms <- if(missing(data)) terms(formula)
           else              terms(formula, data=data)
m$formula <- Terms
m[[1]] <- as.name("model.frame")
m <- eval(m, sys.parent())
mt <- attr(m, "terms")
intercept<-attr(mt, "intercept")
Y <- model.extract(m, "response")
if (!inherits(Y, "Surv")) stop("Response must be a survival object")

if (attr(m[, 1], "type") == "right") {
  time2  <- m[, 1][, "time"]; time   <- rep(0,length(time2));
  status <- m[, 1][, "status"]    } else 
if (attr(m[, 1], "type") == "counting") {
  time   <- m[, 1][,1]; time2  <- m[, 1][,2]; status <- m[, 1][,3]; } else {
  stop("only right-censored or counting processes data") } 

X <- model.matrix(Terms, m)[,-1,drop=FALSE]; covnamesX<-dimnames(X)[[2]]; 
desX<-as.matrix(X);
if(is.matrix(desX) == TRUE) pg <- as.integer(dim(desX)[2])
if(is.matrix(desX) == TRUE) nx <- as.integer(dim(desX)[1])
px<-1; 
Ntimes <- sum(status); 

# adds random noise to make survival times unique
if (sum(duplicated(time2[status==1]))>0) {
#cat("Non unique survival times: break ties ! \n")
#cat("Break ties yourself\n");
ties<-TRUE
dtimes<-time2[status==1]
index<-(1:length(time2))[status==1]
ties<-duplicated(dtimes); nties<-sum(ties); index<-index[ties]
dt<-diff(sort(time2)); dt<-min(dt[dt>0]);
time2[index]<-time2[index]+runif(nties,0,min(0.001,dt/2));
} else ties<-FALSE; 

start<-time; stop<-time2; 
dtimes<-time2[status==1]; 
times<-c(start.time,dtimes[dtimes>start.time]); 
times<-sort(times);
if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
times<-times[times<=maxtimes]
Ntimes <- length(times); 
 
 if ((nrow(X)!=nrow(data)) && (!is.null(id))) stop("Missing values in design matrix not allowed with id\n"); 
### if (nrow(X)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 

########################################################################
if (is.null(id)==TRUE) {antpers<-length(time); id<-0:(antpers-1); }
else { pers<-unique(id); antpers<-length(pers); 
       id<-as.integer(factor(id,labels=1:(antpers)))-1; 
}

cluster.call<-clusters; 
if (is.null(clusters)== TRUE) {clusters<-id; antclust<-antpers;} else {
       clus<-unique(clusters); antclust<-length(clus); 
       clusters <- as.integer(factor(clusters, labels = 1:(antclust))) - 1;
}

if ((!is.null(max.clust))) if (max.clust<antclust) {
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
###	clusters <- as.integer(factor(qqc, labels = 1:max.clust)) -1
	antclust <- max.clust    
  }                


  if (resample.iid == 1) {
    biid <- double(Ntimes* antclust );
    gamiid<- double(antclust *pg);
  } else {
    gamiid <- biid <- NULL;
  }

  if ((length(beta)!=pg) && (is.null(beta)==FALSE)) beta <- rep(beta[1],pg); 
  if ((is.null(beta))) {
        if ( (attr(m[, 1], "type") == "right" ) ) beta<-coxph(Surv(stop,status)~X)$coef
        else beta<-coxph(Surv(start,stop,status)~X)$coef; 
  }

if (residuals==1) {
cumAi<-matrix(0,Ntimes,antpers*1);
cumAiiid<-matrix(0,Ntimes,antpers*1); 
} else { cumAi<-0; cumAiiid<-0; }

cumint<-matrix(0,Ntimes,px+1); 
vcum<-matrix(0,Ntimes,px+1);
Rvcu<-matrix(0,Ntimes,px+1);
score<-beta;
Varbeta<-matrix(0,pg,pg); Iinv<-matrix(0,pg,pg);
RVarbeta<-matrix(0,pg,pg);
if (sim==1) Uit<-matrix(0,Ntimes,50*pg) else Uit<-NULL;

test<-matrix(0,antsim,2*px); testOBS<-rep(0,2*px); unifCI<-c();
testval<-c();
rani<--round(runif(1)*10000); 
Ut<-matrix(0,Ntimes,pg+1); simUt<-matrix(0,antsim,pg);
loglike<-0; 

########################################################################


###cat("Proportional odds model \n"); 
###dyn.load("Gprop-odds.so")


nparout<- .C("transsurv",
as.double(times),as.integer(Ntimes),as.double(desX),
as.integer(nx),as.integer(pg),as.integer(antpers),
as.double(start),as.double(stop), as.double(beta),
as.integer(Nit), as.double(cumint), as.double(vcum),
as.double(Iinv),as.double(Varbeta),as.integer(detail),
as.integer(sim),as.integer(antsim),as.integer(rani),
as.double(Rvcu),as.double(RVarbeta),as.double(test),
as.double(testOBS),as.double(Ut),as.double(simUt),
as.double(Uit),as.integer(id),as.integer(status),
as.integer(weighted.test),as.integer(ratesim),as.double(score),
as.double(cumAi),as.double(cumAiiid),as.integer(residuals),
as.double(loglike),as.integer(profile),as.integer(sym),
as.integer(baselinevar),as.integer(clusters),as.integer(antclust), 
as.double(biid),as.double(gamiid),PACKAGE="timereg");


gamma<-matrix(nparout[[9]],pg,1);
cumint<-matrix(nparout[[11]],Ntimes,px+1);
vcum<-matrix(nparout[[12]],Ntimes,px+1);
Iinv<-matrix(nparout[[13]],pg,pg);
Varbeta<--matrix(nparout[[14]],pg,pg);
Rvcu<-matrix(nparout[[19]],Ntimes,px+1);
RVarbeta<--matrix(nparout[[20]],pg,pg);
score<-matrix(nparout[[30]],pg,1);
Ut<-matrix(nparout[[23]],Ntimes,pg+1);
loglike<-nparout[[34]]

if (residuals==1) {
cumAi<-matrix(nparout[[31]],Ntimes,antpers*1);
cumAiiid<-matrix(nparout[[32]],Ntimes,antpers*1);
cumAi<-list(time=times,dmg=cumAi,dmg.iid=cumAiiid);} else cumAi<-NULL;

 if (resample.iid==1)  {
    biid<-matrix(nparout[[40]],Ntimes,antclust);
    gamiid<-matrix(nparout[[41]],antclust,pg) 
    gamiid  <-  t(Iinv %*% t(gamiid))
    B.iid<-list();
    for (i in (1:antclust)) {
    B.iid[[i]]<-matrix(biid[,i],ncol=1);
    colnames(B.iid[[i]])<-"Baselineiid"; 
    }
      colnames(gamiid)<-covnamesX
  } else B.iid<-gamiid<-NULL;


if (sim==1) {
Uit<-matrix(nparout[[25]],Ntimes,50*pg); UIt<-list();
for (i in (0:49)*pg) UIt[[i/pg+1]]<-as.matrix(Uit[,i+(1:pg)]);
simUt<-matrix(nparout[[24]],antsim,pg);
test<-matrix(nparout[[21]],antsim,2*px); testOBS<-nparout[[22]];
supUtOBS<-apply(abs(as.matrix(Ut[,-1])),2,max);
for (i in 1:(2*px)) testval<-c(testval,pval(test[,i],testOBS[i]));
for (i in 1:px) unifCI<-c(unifCI,percen(test[,i],0.95));
testUt<-c();
for (i in 1:pg) testUt<-c(testUt,pval(simUt[,i],supUtOBS[i]));

pval.testBeq0<-as.vector(testval[1:px]);
pval.testBeqC<-as.vector(testval[(px+1):(2*px)]);
obs.testBeq0<-as.vector(testOBS[1:px]);
obs.testBeqC<-as.vector(testOBS[(px+1):(2*px)]);
sim.testBeq0<-as.matrix(test[,1:px]);
sim.testBeqC<-as.matrix(test[,(px+1):(2*px)]);
sim.supUt<-as.matrix(simUt);
}

if (sim!=1) {
testUt<-NULL;test<-NULL;unifCI<-NULL;supUtOBS<-NULL;UIt<-NULL;testOBS<-NULL;testval<-NULL;
pval.testBeq0<- pval.testBeqC<- obs.testBeq0<- obs.testBeqC<- sim.testBeq0<-
sim.testBeqC<-NULL; testUt<-NULL; sim.supUt<-NULL;
}

ud<-list(cum=cumint,var.cum=vcum,robvar.cum=Rvcu,
gamma=gamma,var.gamma=Varbeta,robvar.gamma=RVarbeta,
resid.dMG=cumAi,D2linv=Iinv,score=score,loglike=loglike,
pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
obs.testBeq0=obs.testBeq0,obs.testBeqC=obs.testBeqC,
sim.testBeq0= sim.testBeq0,sim.testBeqC=sim.testBeqC,
conf.band=unifCI,
test.procProp=Ut,sim.test.procProp=UIt,pval.Prop=testUt,
sim.supProp=sim.supUt,prop.odds=TRUE,gamma.iid=gamiid,B.iid=B.iid)

colnames(ud$cum)<-colnames(ud$var.cum)<- c("time","Baseline")
if (robust==1) colnames(ud$robvar.cum)<- c("time","Baseline"); 

if (px>0) {
if (sim==1) {
colnames(ud$test.procProp)<-c("time",covnamesX)
names(ud$pval.Prop)<-covnamesX
names(ud$conf.band)<-names(ud$pval.testBeq0)<-
names(ud$pval.testBeqC)<-names(ud$obs.testBeq0)<- 
names(ud$obs.testBeqC)<-colnames(ud$sim.testBeq0)<-"Baseline";
} }

rownames(ud$gamma)<-c(covnamesX); colnames(ud$gamma)<-"estimate";
rownames(ud$score)<-c(covnamesX); colnames(ud$score)<-"score";
namematrix(ud$var.gamma,covnamesX);
namematrix(ud$robvar.gamma,covnamesX);
namematrix(ud$D2linv,covnamesX);


attr(ud,"Call")<-call; 
attr(ud,"Formula")<-formula; 
attr(ud,"id")<-id.call; 
attr(ud,"basesim") <- 1
attr(ud,"type") <- "survival"
class(ud)<-"cox.aalen"
return(ud); 
} ## }}} 

