bestfit.search.random<-function(statevar,lagstate,covariate,P,R,Q,indexBGlobal,indexCGlobal,...){
#====================================================================================
# INITIALIZE BEST-FIT MODEL SEARCH (lowest AIC):
#====================================================================================

# Set the random number generator
set.seed(sum(100*proc.time(),na.rm=T))

hitB<-matrix(0,nrow=P,ncol=P)
hitC<-matrix(0,nrow=P,ncol=R)

bestGlobalB<-matrix(0,nrow=P,ncol=P)
bestGlobalC<-matrix(0,nrow=P,ncol=R)

rep.count<-rep(0,P)

cancel<-tclVar(0)
ttinfo<-"-\n-"
ttinfo2<-""
tt<-tktoplevel()
tkwm.geometry(tt,"300x250-30+30")
tkwm.title(tt,"Search Progress")
ttstatus<-tklabel(tt,text=ttinfo)
ttstatus2<-tklabel(tt,text=ttinfo2)
pb<-ttkprogressbar(parent=tt,orient="horizontal",
	length=200,mode="determinate",maximum=P,value=0)
pb2<-ttkprogressbar(parent=tt,orient="horizontal",
	length=100,mode="determinate",maximum=100,value=0)
tkgrid(tklabel(tt,text=""))
tkgrid(ttstatus)
tkgrid(tklabel(tt,text=""))
tkgrid(pb)
tkgrid(pb2)
tkgrid(ttstatus2)
ttmeanreps<-tklabel(tt,text="")
tkgrid(ttmeanreps)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text=" Cancel ",font=tkfont.create(size=10,weight="bold"),
	command=function() {tclvalue(cancel)<-1}))
tkgrid.columnconfigure(tt,pb,minsize=300)

all.models<-vector("list",P)

for(dv in 1:P){	#  *************** START Main loop ***************
# (for each of the variates...)

ttinfo<-paste("Searching for best-fit model for\nvariate",dv,"of",P,"...")
tkconfigure(pb,value=(dv-1))
tkconfigure(ttstatus,text=ttinfo)

indexB<-indexBGlobal[dv,]
indexC<-indexCGlobal[dv,]

bestAIC<-1000000
bestBIC<-bestAIC

# reset counters to eliminate coefficients with < 15 hits
countB<-c(P,0)
countC<-c(R,0)

while((countB[1]+countC[1])!=(countB[2]+countC[2])){	##### < 15 hit loop #####

# count the number of times models are constructed
rep.count[dv]<-rep.count[dv]+1

ttinfo2<-paste("search iteration: ",rep.count[dv])
tkconfigure(ttstatus2,text=ttinfo2)

hitB[dv,]<-0
hitC[dv,]<-0

for(globalindex in 1:100){	# ^^^^^^ START 100 Best loop ^^^^^^
# (find bestAIC of 100 lowestAIC models)

lowestAIC<-1000000
lowestBIC<-lowestAIC

for(rept in 1:100){			# ------ START 100 Random loop ------
# (find lowestAIC of 100 random models)

E<-matrix(0,nrow=Q,ncol=1)	# residuals
A<-matrix(0,nrow=1,ncol=1)	# intercepts for the variates
B<-matrix(0,nrow=1,ncol=P)	# parameters for the variates
C<-matrix(0,nrow=1,ncol=R)	# parameters for the covariates
Yhat<-matrix(0,nrow=Q,ncol=1)	# estimates

# varind = matrix for the variates : which interactions to include
vareleB<-which(indexB==.5)
randvec<-round(runif(length(vareleB),0,1))
varind<-indexB
varind[varind==.5]<-randvec

# covarind = matrix for the covariates : which interactions to include
vareleC<-which(indexC==.5)
randvec<-round(runif(length(vareleC),0,1))
covarind<-indexC
covarind[covarind==.5]<-randvec

Y<-statevar[,dv]						# dependent variate
lv<-which(varind==1)
cv<-which(covarind==1)

nvar<-sum(varind)

X<-cbind(1,lagstate[,lv],covariate[,cv])		# predictors

# least squares estimates
beta<-solve( (t(X)%*%X),(t(X)%*%Y) )

# calculate residuals for the dependent variate
Yhat<-X%*%beta
E[,1]<-Y-Yhat

# parameter estimates
A<-beta[1,1]
if (length(lv)>0) B[1,lv]<-beta[2:(nvar+1),1]
if (length(cv)>0) C[1,cv]<-beta[(nvar+2):length(beta[,1]),1]

# calculate log-likelihood
sigma<-t(E)%*%E/Q
lnlike<- -Q*(1/2)*log(2*pi)-(Q/2)*log(det(sigma))-Q/2

par<-sum(A!=0,B!=0,C!=0)+1
SSAIC<- -2*lnlike+2*par
SSBIC<- -2*lnlike+par*log(Q)

if(SSAIC<lowestAIC){
	lowestAIC<-SSAIC
	lowestB<-B
	lowestC<-C}

}						# ------- END 100 Random loop -------

all.models[[dv]]<-rbind(all.models[[dv]],c(rep.count[dv],lowestAIC,lowestB,lowestC))

if(lowestAIC<bestAIC){
	bestAIC<-lowestAIC
	bestB<-lowestB
	bestC<-lowestC}

patternB<-lowestB!=0
patternC<-lowestC!=0
hitB[dv,]<-hitB[dv,]+patternB
hitC[dv,]<-hitC[dv,]+patternC

tkconfigure(pb2,value=globalindex)
if(as.logical(as.numeric(tclvalue(cancel)))) break()
}					# ^^^^^^^ END 100 Best loop ^^^^^^^

indexB<-indexB*(hitB[dv,]>=15)
indexC<-indexC*(hitC[dv,]>=15)

countB<-c(sum(indexB!=0),countB)
countC<-c(sum(indexC!=0),countC)

if(as.logical(as.numeric(tclvalue(cancel)))) break()
}							##### END < 15 hit loop #####

meanreps<-round(mean(rep.count[rep.count!=0]),2)
tkconfigure(ttmeanreps,
	text=paste("( mean: ",meanreps,")"))

bestGlobalB[dv,]<-bestB
bestGlobalC[dv,]<-bestC

if(as.logical(as.numeric(tclvalue(cancel)))) break()
}			#  **************** END Main loop ****************

if(as.logical(as.numeric(tclvalue(cancel)))) {tkdestroy(tt);NULL} else{
tkconfigure(pb,value=dv)
ttinfo<-"COMPLETE\n - "
ttinfo2<-"-"
tkconfigure(ttstatus,text=ttinfo)
tkconfigure(ttstatus2,text=ttinfo2)
tkconfigure(ttmeanreps,text=ttinfo2)
Sys.sleep(2)
tkdestroy(tt)

for(i in 1:length(all.models)) all.models[[i]]<-all.models[[i]][,-1]

c(
list(bestGlobalB=bestGlobalB,bestGlobalC=bestGlobalC),
all.models=list(all.models)
)
}
}