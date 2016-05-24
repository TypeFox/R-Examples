pssm.survivalcurv <-
function(x,cov1,cov2,timeToProgression=FALSE,covariance=TRUE){
#library(abind)
#input is a pssm object and output is the survival curve or time to progression function
#outputs a function that takes a vector t of length L and outputs an 
#M X L matrix of values with M being the number of rows of COV1,and COV2,
#if covariance=TRUE the output of the function is a list consisting of the estimates and it's variance covariance matrix.
#curv can be survival,progression(for survivorship) and dsurvival,dprogression for density. 
#the outputed function takes on an optional vector arguement of the prior 
#precision of the parameters  
# The parameters are
# first the log-hazards, then progression coefficients and then survival coefficients

sp=length(x@formula.progression)==0
ss=length(x@formula.survival)==0
cd1=(!is.null(cov1))&(!sp)
cd2=(!is.null(cov2))&(!ss)
if(cd1) dc1=dim(cov1) else dc1=rep(0,2)
if(cd2) dc2=dim(cov2) else dc2=rep(0,2)
#haz=matrix(exp(x@hazard.progression),1,x@intervals)
#chaz=matrix(cumsum(haz),dc1[1],x@intervals,byrow=TRUE)
namesc1=c("rep","tdeath","cdeath","tprog0","tprog1",x@progression.covariate.list,x@survival.covariate.list)
namesc=namesc1[namesc1!=""]

if(!(sp|ss)){
	#tfcn1<-llikef(x@progression.covariate.list,x@survival.covariate.list,dt,m=x@intervals,accumulate=FALSE,gradient=TRUE)
  tfcn2<-function(dt) llikef(x@progression.covariate.list,x@survival.covariate.list,dt,
                             m=x@intervals,accumulate=FALSE,gradient=FALSE)
  tfcn1<-function(dt) {
    lp=dim(dt)[1]
    
    outp<-function(vv){
      out1<-tfcn2(dt)(vv)
	    gradient=matrix(NaN,lp,length(vv))
	   for(i in (1:lp)){
	     try(gradient[i,]<-grad(function (xt) tfcn2(dt[i,])(xt),vv),silent=FALSE)    
       }
     attr(out1,"gradient")<-gradient
 
	   return(out1)
     }
  return(outp)
  }

}	else {
	if (sp){ 
	tfcn1<-function(dt) {
	
		lp=dim(dt)[1]
		outp=function(vv){
		out1<-rsurv(x@survival.covariate.list,dt,m=x@intervals,accumulate=FALSE)(vv)
		gradient=matrix(NaN,lp,length(vv))
		for (i in (1:lp)) 
		try(gradient[i,]<-grad(function(xt) rsurv(x@survival.covariate.list,data.frame(dt[i,]),m=x@intervals,accumulate=FALSE)(xt),vv),silent=TRUE)
		attr(out1,"gradient")<-gradient

    return(out1)
	   }
      return(outp)
	}
    	tfcn2<-function(dt) rsurv(x@survival.covariate.list,dt,m=x@intervals,accumulate=FALSE)
} else
    {
      tfcn1<-function(dt){
    
	   mtt=dim(dt)[1]
	   outp=function(vv){

	   out1=rprog(x@progression.covariate.list,dt,m=x@intervals,accumulate=FALSE)(vv)
	   gradient=matrix(NaN,mtt,length(vv))
       for (i in (1:mtt))
	   try(gradient[i,]<-grad( function(xt) rprog(x@progression.covariate.list,data.frame(dt[i,]),
                                                m=x@intervals,accumulate=FALSE)(xt),vv),silent=TRUE)
	   attr(out1,"gradient")<-gradient
	   return(out1)}
	    return(outp)
		}
    tfcn2<-function(dt) rprog(x@progression.covariate.list,dt,m=x@intervals,accumulate=FALSE)
      }
}

#this is the workhorse function for calling llikef

#table of calls
curv1=NULL
#expands the list of curves to draw
if (!(sp|ss)) curv1=c('s1','s2') else curv1='s2'
#for (i in 1:length(curv)) if (sp) curv1=c(curv1,c('s1','s2')) else curv1=c(curv1,substr(curv[i],1,2))
curv1=unique(curv1)
#forms data frame with data
nt=function(t1,t2,t3,t4,t5) return(data.frame(rep=t1,tdeath=t2,cdeath=t3,tprog0=t4,tprog1=t5,stringsAsFactors=FALSE))
#repeats matrix m times
repframe=function(x,m){out=NULL
	for (i in (1:m)) out=rbind(out,x)
return(out)}

#used with nt to create data frames
#pr,#s2 probability that patient had progression after time t
#s1 probablity that they had pregression before t but lived after t
#dp density for progression at time t, given alive after t
#ds density of death at t, with progression between 0 and t
cls=function(t){
	ls=list(pr=nt('pr',t,0,t,NA),s1=nt('s1',t,0,0,t),s2=nt('s2',t,0,t,NA),dp=nt('dp',t,0,t,t),ds=nt('ds',t,1,0,t))
	#currently only uses s1 and s2
	inp=list(NULL,NULL)
	lt=length(t)
	for (i in 1:length(curv1)) {
		inp[[1]]=rbind(inp[[1]],ls[[curv1[i]]])
        inp[[2]]=rbind(inp[[2]],diag(rep(1,lt)))}
    return(inp)
}
idt=function(ts){
   tt=ts[[1]]
	vs=dim(tt)[1]
	rv=max(max(dc1[1],dc2[1]),1)
	mm=rv*vs

	inp=data.frame(rep=rep(tt[,1],rv),time=rep(tt[,2],rv),stringsAsFactors=FALSE)
	if (cd1){		fcov1=data.frame(matrix(rep(cov1,each=vs),mm,dc1[2]));	names(fcov1)<-x@progression.covariate.list
		inp=cbind(inp,fcov1)}
	if (cd2){		fcov2=data.frame(matrix(rep(cov2,each=vs),mm,dc2[2]));	names(fcov2)<-x@survival.covariate.list
		inp=cbind(inp,fcov2)}
	rownames(inp)<-as.character(1:mm)
	inp2=diag(rep(1,rv))%x%ts[[2]]
	return(list(inp,inp2))
}
sortrows<-function(x,r){
	#sort data frame on first m rows
	out=x
	for (i in seq(r,1,-1)) out=out[order(out[,i]),]
	return(out)
}
outp=function(ts,p){
	#tt vs+1 x 4k matrix depending on the function tdeath,cdeath,tprog0,tprog1
	#output is number of covariates groups X (VS+ number of covariate

	tt=ts[[1]]
	vs=dim(tt)[1]  #number of timepoints
	rv=max(max(dc1[1],dc2[1]),1)      #number of curves #dc1 #of covariates
	mm=rv*vs   #size of data frame needed
dt=data.frame(repframe(tt,rv),idt(ts)[[1]][,-(1:2)])
	#dt=data.frame(cbind(rep(tt[,1],mm),rep(tt[,2],mm),rep(tt[,3],mm),rep(tt[,4],mm),rep(tt[,5],mm),idt(tt)[,3:(dc1[2]+dc2[2]+2)]))
    #dt=data.frame(tt,inp[,c(x@progression.covariate.list,x@survival.covariate.list)])
    names(dt)<-namesc
 	if(covariance){
    fcc=tfcn1(dt)

		tmp1=fcc(x@estimates)
		fc=exp(tmp1)
    me=length(x@estimates)
    dlogf=attr(tmp1,"gradient")
		nrows=dim(dlogf)[1]
    vlogf=dlogf%*%x@covariance.estimates%*%t(dlogf)
  
    if (p!=0){
      #new way 
      pp=p
      E12=matrix(x@covariance.estimates[1:(me-1),me],me-1,1)
      e22=x@covariance.estimates[me,me]
      ED=rbind(cbind(matrix(0,me-1,me-1),pp*E12),
               cbind(matrix(0,1,me-1),pp*e22))
      EDE=rbind(cbind(pp*E12%*%t(E12),e22*pp*E12),
                cbind(e22*pp*t(E12),pp*e22^2))
      pest=x@estimates-(ED%*%matrix(x@estimates,me,1))/(1+e22*pp)
      cvar=x@covariance.estimates-EDE/(1+pp*e22)
      tmp1=fcc(pest)
      fc=exp(tmp1)
      dlogf=attr(tmp1,"gradient")
      vlogf=dlogf%*%cvar%*%t(dlogf)
     
       }   
		  fk=matrix(fc,mm,mm)
		  cdc=fk*vlogf*t(fk)

ou=cbind(matrix(fc,nrows,1),cdc)
colnames(ou)<-c('estimate',as.character(1:mm))
rownames=as.character(1:mm)
return(ou)    
	}
		else {
   
		  fcc=tfcn2(dt)
	    fc=exp(fcc(x@estimates))
		  ou=matrix(fc,length(fc),1)
		  colnames(ou)<-c('estimate')
		  rownames(ou)<as.character(1:length(fc))
		  return(ou)}
}
out<-function(ts,p=0) {
#p=precision of baysean
#patch for t=0;
rt<-function(t){ifelse(t==0,0.00001,t*x@rescale)}
rti<-function(t){ifelse(t==0.00001,0,t/x@rescale)}
t=rt(ts)
vin=cls(t)
vt=idt(vin)
vu=outp(vin,p)
xvals=colnames((vt)[[1]])[-1]
vt[[1]][,2]<-rti(vt[[1]][,2])
v=cbind(vt[[1]],vu)
ninfo=dc1[2]+dc2[2]+2
names(v)<-c(names(vt[[1]]),colnames(vu))
if (!(sp|ss|timeToProgression)){
#v=sortrows(v,ninfo)
#vs=merge(v[v$rep=="s1",],v[v$rep=="s2",],by=xvals)
#nhalf=dim(v)[1]/2
#m1=diag(rep(1,nhalf))
#m1=rbind(m1,m1)
#ou=cbind(v[1:nhalf,1:ninfo],matrix(matrix(v[,"estimate"],1,2*nhalf)%*%m1,nhalf,1))
#needs to be fixed at some point if we do more than survival
ou=cbind(v[v[,1]==v[1,1],1:ninfo],estimate=t(vt[[2]])%*%matrix(v[,"estimate"],length(v[,"estimate"]),1))
names(ou)[[ninfo+1]]="estimate"
}
#ou=data.frame(vs[,xvals,drop=FALSE],estimate=vs[,"estimate.x"]+vs[,"estimate.y"],stringsAsFactors=FALSE)}
else {
ou=data.frame(v[,c("rep",xvals),drop=FALSE],estimate=v[,"estimate"],stringsAsFactors=FALSE)
}	

if (covariance){
	cov=t(vt[[2]])%*% as.matrix(v[,-(1:(ninfo+1))]) %*%vt[[2]]
	attr(ou,"covariance")<-cov
	}
return(ou)}
return(out)
}