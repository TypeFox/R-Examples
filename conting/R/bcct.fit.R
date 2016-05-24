bcct.fit <-
function(priornum,maximal.mod,IP,eta.hat,ini.index,ini.beta,ini.sig,iters,save,name,null.move.prob,a,b,progress){

if(is.null(name)){
name_RJACC<-"RJACC.txt"
name_MHACC<-"MHACC.txt"
name_BETA<-"BETA.txt"
name_MODEL<-"MODEL.txt"
name_SIG<-"SIG.txt"} else{
name_RJACC<-paste(name,"RJACC.txt",sep="")
name_MHACC<-paste(name,"MHACC.txt",sep="")
name_BETA<-paste(name,"BETA.txt",sep="")
name_MODEL<-paste(name,"MODEL.txt",sep="")
name_SIG<-paste(name,"SIG.txt",sep="")}       					## names for files if we are saving them

big.X<-maximal.mod$x								## maximal design matrix			
y<-maximal.mod$y									## responses
data<-maximal.mod$data								## data

curr.index<-ini.index								## current index (binary vector)
curr.X<-big.X[,curr.index==1]							## current design matrix
curr.ivar<-IP[curr.index==1,curr.index==1]					## current prior scale matrix
MODEL<-c()					## current model (hexidecimal)
BETA<-c()				## matrix for betas
curr.beta<-ini.beta								## current betas
SIG<-c()									## current prior variance
curr.sig<-ini.sig								## ditto
rj_acc<-c()									## setting up vector for RJ accept/rejects	
mh_acc<-c()									## setting up vector for MH accept/rejects
counter<-0									## iteration counter
if(progress){									
pb<-txtProgressBar(min = 0, max = iters, style = 3)}		## set up progress bar
while(counter<iters){

prop_a_model<-prop_mod(curr.index=curr.index,data=data,maximal.mod=maximal.mod,null.move.prob=null.move.prob)	## propose a model
prop.index<-prop_a_model$new.index						## proposal index

if(prop_a_model$type=="null"){									## stay in same model
new.beta<-iwls_mh(curr.y=y,curr.X=curr.X,curr.beta=curr.beta,iprior.var=curr.ivar/curr.sig)	## do MH
if(all(new.beta==curr.beta)){mh_acc<-c(mh_acc,0)} else{mh_acc<-c(mh_acc,1)} 			## did we move?
new.index<-curr.index}										## new index (same as last)

if(prop_a_model$type!="null"){									## proposing new model

rho_bot<-(1-prop_a_model$null.move.prob)/prop_a_model$total.choices				## proposal prob (current to proposed)
prop_a_rev<-prop_mod(curr.index=prop.index,data=data,maximal.mod=maximal.mod,null.move.prob=0)				## proposal prob (proposed to current)	
rho_top<-(1-prop_a_model$null.move.prob)/prop_a_rev$total.choices					## ditto

rj<-RJ_update(prop.index=prop.index,curr.index=curr.index,curr.beta=curr.beta,eta.hat=eta.hat,
curr.y=y,big.X=big.X,proposal.probs=c(rho_top,rho_bot),i.prop.prior.var=IP[prop.index==1,prop.index==1]/curr.sig,
i.curr.prior.var=curr.ivar/curr.sig)								## do RJ

new.beta<-rj$new.beta
new.index<-rj$new.index										## new betas and index
if(all(curr.index==new.index)){rj_acc<-c(rj_acc,0)} else{rj_acc<-c(rj_acc,1)}}			## did we move?

if(priornum==2){										## update sig if we have an SBH prior
#a<-0.001
#b<-0.001
iR<-IP[new.index==1,new.index==1]
curr.sig<-1/rgamma(n=1,shape=0.5*(length(new.beta)-1+a),rate=0.5*(b+as.vector(matrix(new.beta[-1],nrow=1)%*%iR[-1,-1]%*%matrix(new.beta[-1],ncol=1))))
}

curry<-rep(0,dim(big.X)[2])
curry[new.index==1]<-new.beta

curr.index<-new.index
curr.beta<-new.beta
curr.X<-big.X[,curr.index==1]
curr.ivar<-IP[curr.index==1,curr.index==1]
BETA<-rbind(BETA,curry)
SIG<-c(SIG,curr.sig)
MODEL<-c(MODEL,index2model(new.index))								## update new parameters

counter<-counter+1
if(progress){
setTxtProgressBar(pb, counter)}									## update progress bar

if(save>0){
if(counter%%save==0){

write.table(file=name_BETA,x=BETA,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_MODEL,x=MODEL,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_SIG,x=SIG,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_RJACC,x=rj_acc,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_MHACC,x=mh_acc,row.names=FALSE,col.names=FALSE,append=TRUE)		## save progress

rj_acc<-c()
mh_acc<-c()
BETA<-c()
MODEL<-c()
SIG<-c()}}}
if(progress){
close(pb)}													## close progress bar

list(BETA=BETA,SIG=SIG,MODEL=MODEL,rj_acc=rj_acc,mh_acc=mh_acc)}
