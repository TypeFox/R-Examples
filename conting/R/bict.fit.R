bict.fit <-
function(priornum,missing1,missing2,maximal.mod,IP,eta.hat,ini.index,ini.beta,ini.sig,ini.y0,iters,save,name,null.move.prob,a,b,progress){

missing<-c(missing1,missing2)

if(is.null(name)){
name_RJACC<-"RJACC.txt"
name_MHACC<-"MHACC.txt"
name_BETA<-"BETA.txt"
name_MODEL<-"MODEL.txt"
name_SIG<-"SIG.txt"
name_Y0<-"Y0.txt"} else{
name_RJACC<-paste(name,"RJACC.txt",sep="")
name_MHACC<-paste(name,"MHACC.txt",sep="")
name_BETA<-paste(name,"BETA.txt",sep="")
name_MODEL<-paste(name,"MODEL.txt",sep="")
name_SIG<-paste(name,"SIG.txt",sep="")
name_Y0<-paste(name,"Y0.txt",sep="")}

big.X<-maximal.mod$x
y<-maximal.mod$y
data<-maximal.mod$data

curr.index<-ini.index
curr.X<-big.X[,curr.index==1]
curr.ivar<-IP[curr.index==1,curr.index==1]
#MODEL<-as.character(index2model(curr.index))
MODEL<-c()
#BETA<-matrix(0,nrow=1,ncol=dim(big.X)[2])
#BETA[,curr.index==1]<-ini.beta
BETA<-c()
curr.beta<-ini.beta
#SIG<-ini.sig
SIG<-c()
curr.sig<-ini.sig
#Y0<-matrix(ini.y0,nrow=1)
Y0<-c()
curr.y0<-ini.y0
rj_acc<-c()
mh_acc<-c()
counter<-0
if(progress){									
pb<-txtProgressBar(min = 0, max = iters, style = 3)}		## set up progress bar
while(counter<iters){

curr.y<-y
curr.y[missing]<-curr.y0

prop_a_model<-prop_mod(curr.index=curr.index,data=data,maximal.mod=maximal.mod,null.move.prob=null.move.prob)
prop.index<-prop_a_model$new.index

if(prop_a_model$type=="null"){
new.beta<-iwls_mh(curr.y=curr.y,curr.X=curr.X,curr.beta=curr.beta,iprior.var=curr.ivar/curr.sig)
if(all(new.beta==curr.beta)){mh_acc<-c(mh_acc,0)} else{mh_acc<-c(mh_acc,1)} 
new.index<-curr.index}

if(prop_a_model$type!="null"){

rho_bot<-(1-prop_a_model$null.move.prob)/prop_a_model$total.choices
prop_a_rev<-prop_mod(curr.index=prop.index,data=data,maximal.mod=maximal.mod,null.move.prob=0)
rho_top<-(1-prop_a_model$null.move.prob)/prop_a_rev$total.choices

rj<-RJ_update(prop.index=prop.index,curr.index=curr.index,curr.beta=curr.beta,eta.hat=eta.hat,
curr.y=curr.y,big.X=big.X,proposal.probs=c(rho_top,rho_bot),i.prop.prior.var=IP[prop.index==1,prop.index==1]/curr.sig,
i.curr.prior.var=curr.ivar/curr.sig)

new.beta<-rj$new.beta
new.index<-rj$new.index
if(all(curr.index==new.index)){rj_acc<-c(rj_acc,0)} else{rj_acc<-c(rj_acc,1)}}

if(priornum==2){
#a<-0.001
#b<-0.001
iR<-IP[new.index==1,new.index==1]
curr.sig<-1/rgamma(n=1,shape=0.5*(length(new.beta)-1+a),rate=0.5*(b+as.vector(matrix(new.beta[-1],nrow=1)%*%iR[-1,-1]%*%matrix(new.beta[-1],ncol=1))))
}

mutar<-exp(as.vector(big.X[,new.index==1]%*%matrix(new.beta,ncol=1)))
curr.y0<-rep(0,length(missing))
curr.y0[1:length(missing1)]<-rpois(n=length(missing1),lambda=mutar[missing1])
if(length(missing2)>0){
for(i in 1:length(missing2)){
ppp<-dpois(x=1:y[missing2[i]],lambda=mutar[missing2[i]],log=TRUE)
ppp<-ppp-max(ppp)
curr.y0[length(missing1)+i]<-sample(x=1:y[missing2[i]],size=1,prob=exp(ppp))}}

curry<-rep(0,dim(big.X)[2])
curry[new.index==1]<-new.beta

curr.index<-new.index
curr.beta<-new.beta
curr.X<-big.X[,curr.index==1]
curr.ivar<-IP[curr.index==1,curr.index==1]
BETA<-rbind(BETA,curry)
SIG<-c(SIG,curr.sig)
MODEL<-c(MODEL,index2model(new.index))
Y0<-rbind(Y0,curr.y0)

counter<-counter+1
if(progress){
setTxtProgressBar(pb, counter)}									## update progress bar

if(save>0){
if(counter%%save==0){

write.table(file=name_BETA,x=BETA,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_MODEL,x=MODEL,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_SIG,x=SIG,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_Y0,x=Y0,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_RJACC,x=rj_acc,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_MHACC,x=mh_acc,row.names=FALSE,col.names=FALSE,append=TRUE)

rj_acc<-c()
mh_acc<-c()
BETA<-c()
MODEL<-c()
SIG<-c()
Y0<-c()}}}
if(progress){
close(pb)}

list(BETA=BETA,SIG=SIG,MODEL=MODEL,Y0=Y0,rj_acc=rj_acc,mh_acc=mh_acc)}
