pssm.simulate <-
function(nruns,theta1=.2,theta2=.2,phaz.progression=log(-log(.3)/4)*rep(1,5),
	phaz.survival=log(-log(.15)/4)*rep(1,15),accrual=3,followup=2,
  m=5,n=500,rescale=1,etime=4.5,seed=NULL){
#vs is full model, #vsp is model with progression only #vst is survival model only

ests=matrix(NaN,nruns,12)
se.ests=matrix(NaN,nruns,12)
namesc=c("ProgParam","SurvParm",paste("DifferenceAt",etime,sep="_"),
	     "ProgParam.NoSurv","NA",paste("DifferenceAt",etime,sep="_"),
		 "NA","TumorFreeDiff.Param",paste("DifferenceAt",etime,sep="_"),
		 paste("RawDifferenceAt",etime,sep="_"),"ProportionDeaths","ProportionProgression")	
colnames(ests)<-namesc
colnames(se.ests)<-namesc
vs=list()
tm=matrix(c(1,-1),1,2)
stt=0
ust=quote({
	if (!is.null(vtt)){
	if(vtt@convergence==0){
	st=stt*3
	r1=st+(1:2)
    ests[i,r1]=c(vtt@estimates.progression,vtt@estimates.survival)
    se.ests[i,r1]=c(vtt@se.estimates.progression,vtt@se.estimates.survival)
    cc1=pssm.survivalcurv(vtt,cov1=matrix(c(0,1),2,1),cov2=matrix(c(0,1),2,1),covariance=TRUE)(etime)
	r2=st+3
	ests[i,r2]=cc1$estimate[1]-cc1$estimate[2]
	se.ests[i,r2]=sqrt(tm%*%attr(cc1,"covariance")%*%t(tm))
	vs[[i]][1+stt]=vtt
}}}
)
	
for (i in (1:nruns)){
	vs[[i]]=list()
	u=pssm.generate.data(theta1,theta2,phaz.progression,phaz.survival,accrual,followup,m,n,seed=seed)
	#analyses using the full model 
    vtt<-pssm(surv(tprog0,tprog1)~rx,surv(tdeath,cdeath)~rx,u,intervals=m,rescale=rescale)
    stt=0
	eval(ust)
	#analysis of progression only
	vtt=pssm(surv(tprog0,tprog1)~rx,NULL,u,intervals=m,rescale=rescale)
	stt=1
	eval(ust)
	#analysis of tumor free progression
	cdeath1=as.numeric(!is.na(u$tprog1))
    tdeath1=ifelse(is.na(u$tprog1),u$tprog0,u$tprog1)
	vtt=pssm(NULL,surv(tdeath1,cdeath1)~rx,data.frame(tdeath1,cdeath1,rx=u$rx),intervals=m,rescale=rescale)
	stt=2
	eval(ust)
    #analysis of survival	
	vtt0=pssm(NULL,surv(tdeath,cdeath)~1,u[u$rx==0,],intervals=m,rescale=rescale)
	cc0=pssm.survivalcurv(vtt0,NULL,NULL,covariance=TRUE)(etime)
	vtt1=pssm(NULL,surv(tdeath,cdeath)~1,u[u$rx==1,],intervals=m,rescale=rescale)
	cc1=pssm.survivalcurv(vtt1,NULL,NULL,covariance=TRUE)(etime)
	ests[i,10]=cc0$estimate[1]-cc1$estimate[1]
    se.ests[i,10]=sqrt(attr(cc0,"covariance")+attr(cc1,"covariance"))
	vs[[i]][4]=vtt0
	vs[[i]][5]=vtt1
	ests[i,11]=mean(u$cdeath)
	ests[i,12]=mean(!is.na(u$tprog1))
}	
return(list(objects=vs,ests=ests,se.ests=se.ests) )
}
