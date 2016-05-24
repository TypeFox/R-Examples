pois.glm.mix <-
function(reference,response,L,m,max.iter,Kmin,Kmax,m1,m2,t1,t2, msplit,tsplit,mnr){


x<-reference
y<-response
maxnr = mnr
if (is.vector(x)==T) x<-array(x,dim=c(length(x),1))
if (is.vector(y)==T) y<-array(y,dim=c(length(y),1))
if (is.data.frame(x)==T){x<-as.matrix(x);colnames(x)=NULL}
if (is.data.frame(y)==T){y<-as.matrix(y);colnames(y)=NULL}

q<-length(L)
tau<-dim(x)[2]

# control if an extra R graphics device should be opened for printing the clusters during the execution of the algorithm
dn<-T
if (tau>1){
dn<-F
}
if(q>5){
dn<-F
}

dev.new()

if (dn==T){
dev.new()
}

# this is for R CMD check warnings
prev.z<-prev.clust<-prev.alpha<-prev.beta<-0



bbs<-array(data = NA, dim = c(Kmax,3))
runs<-vector("list",length = Kmax)

nc<-Kmin
m=floor(m)
if (m < 1 || m>3) stop("m should be in the set {1,2,3}")
index<-numeric(q)
index[1]<-1
if (q > 1){
index[2:q] <- 1+cumsum(L[1:(q-1)])
}

#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
if (m==1){     #
print(paste("fitting the b_jk parameterization"))#
print(paste("***************************************"))#
print(paste("*                K =", nc,"               *"))#
print(paste("***************************************"))#
run<-bjkmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z,prev.clust,start.type=1)#
icl<-run$icl#
bbs[nc,]<-c(run$bic,run$icl,run$ll)#
runs[[nc]]<-run#
for (nc in (Kmin+1):Kmax){#
print(paste("***************************************"))#
print(paste("*                K =", nc,"               *"))#
print(paste("***************************************"))#
z<-run$z#
ml<-length(run$psim)/(nc-1)#
alpha<-array(run$alpha[ml,,],dim=c(q,nc-1))#
beta<-array(run$beta[ml,,,],dim=c(q,nc-1,tau))#
clust<-run$clust#
run <- tryCatch(
{pissa <- bjkmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z=z,prev.clust=clust,start.type=2,prev.alpha=alpha,prev.beta=beta)},
error = function(pissa){
	pissa<-bjkmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z=z,prev.clust=clust,start.type=1,prev.alpha=alpha,prev.beta=beta);return(pissa)
},
finally = pissa
)
bbs[nc,]<-c(run$bic,run$icl,run$ll)#
runs[[nc]]<-run#
plot(c(1,nc),c(min(bbs[1:nc,1],na.rm=T),max(bbs[1:nc,2],na.rm=T)),type="n",xlab="K",ylab="criterion")#
lines(bbs[1:nc,1])#
lines(bbs[1:nc,2],col=2)#
legend("topright",col=c(1,2),lty=1,c("BIC","ICL"))#
#
# graphics#
if (dn==T){#
if(run$icl<icl){#
icl<-run$icl#
print(paste("found new best model: K=",nc))#
dev.set()#
#############################################################
yy<-array(data = NA, dim = c(dim(y)[1],q))#
for (j in 1:q)yy[,j]<-rowSums(array(y[,index[j]:(index[j]+L[j]-1)],dim=c(dim(y)[1],L[j])))#
gamma<-run$gamma#
ml<-length(run$psim)/nc#
alpha<-array(run$alpha[ml,,],dim=c(q,nc))#
beta<-array(run$beta[ml,,,],dim=c(q,nc))#
op<-par(mfrow=c(q,2))#
clust<-run$clust#
for(j in 1:q){#
yyy<-log(yy[,j])#
plot(as.data.frame(cbind(x,exp(yyy))),col=run$clust,ylab=paste("sum(", L[j],"replicates)"),xlab="x",main=paste("condition: ", j))#
for (k in 1:nc){#
norm <- function(u) exp(alpha[j,k] + beta[j,k]*u + log(sum(exp(gamma[j,1:L[j]]))))#
nx <- length(x[clust==k,1])
if (nx>0){
mm<-max(x[clust==k,1])#
curve(norm,from = 0, to = mm+0.5,col = k, lty=k,lwd=1, add = T)}}#
#
plot(as.data.frame(cbind(x,yyy)),col=run$clust,ylab=paste("log(sum(", L[j],"replicates))"),xlab="x",main=paste("condition: ", j))#
for (k in 1:nc){#
norm <- function(u) alpha[j,k] + beta[j,k]*u + log(sum(exp(gamma[j,1:L[j]])))#
nx <- length(x[clust==k,1])
if (nx>0){
mm<-max(x[clust==k,1])#
curve(norm,from = 0, to = mm+0.5,col = k, lty=k,lwd=1, add = T)}}#
}#
dev.set()#
}#
#############################################################
}#
}#
}#
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################




#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
if (m==2){#
print(paste("fitting the b_j parameterization"))#
print(paste("***************************************"))#
print(paste("*                K =", nc,"               *"))#
print(paste("***************************************"))#
run<-bjmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z,prev.clust,start.type=1)#
icl<-run$icl#
bbs[nc,]<-c(run$bic,run$icl,run$ll)#
runs[[nc]]<-run#
for (nc in (Kmin+1):Kmax){#
print(paste("***************************************"))#
print(paste("*                K =", nc,"               *"))#
print(paste("***************************************"))#
z<-run$z#
ml<-length(run$psim)/(nc-1)#
alpha<-array(run$alpha[ml,,],dim=c(q,nc-1))#
beta<-array(run$beta[ml,,],dim=c(q,tau))#
clust<-run$clust#
run <- tryCatch(
{pissa <- bjmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z=z,prev.clust=clust,start.type=2,prev.alpha=alpha,prev.beta=beta)},
error = function(pissa){
	pissa<-bjmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z=z,prev.clust=clust,start.type=1,prev.alpha=alpha,prev.beta=beta);return(pissa)
},
finally = pissa
)
bbs[nc,]<-c(run$bic,run$icl,run$ll)#
runs[[nc]]<-run#
plot(c(1,nc),c(min(bbs[1:nc,1],na.rm=T),max(bbs[1:nc,2],na.rm=T)),type="n",xlab="K",ylab="criterion")#
lines(bbs[1:nc,1])#
lines(bbs[1:nc,2],col=2)#
legend("topright",col=c(1,2),lty=1,c("BIC","ICL"))#
#
# graphics#
if (dn==T){#
if(run$icl<icl){#
icl<-run$icl#
print(paste("found new best model: K=",nc))#
dev.set()#
#############################################################
yy<-array(data = NA, dim = c(dim(y)[1],q))#
for (j in 1:q)yy[,j]<-rowSums(array(y[,index[j]:(index[j]+L[j]-1)],dim=c(dim(y)[1],L[j])))#
gamma<-run$gamma#
ml<-length(run$psim)/nc#
alpha<-array(run$alpha[ml,,],dim=c(q,nc))#
beta<-array(run$beta[ml,,],dim=c(q,1))#
op<-par(mfrow=c(q,2))#
clust<-run$clust#
for(j in 1:q){#
yyy<-log(yy[,j])#
plot(as.data.frame(cbind(x,exp(yyy))),col=run$clust,ylab=paste("sum(", L[j],"replicates)"),xlab="x",main=paste("condition: ", j))#
for (k in 1:nc){#
norm <- function(u) exp(alpha[j,k] + beta[j,1]*u + log(sum(exp(gamma[j,1:L[j]]))))#
nx <- length(x[clust==k,1])
if (nx>0){
mm<-max(x[clust==k,1])#
curve(norm,from = 0, to = mm+0.5,col = k, lty=k,lwd=1, add = T)}}#
#
plot(as.data.frame(cbind(x,yyy)),col=run$clust,ylab=paste("log(sum(", L[j],"replicates))"),xlab="x",main=paste("condition: ", j))#
for (k in 1:nc){#
norm <- function(u) alpha[j,k] + beta[j,1]*u + log(sum(exp(gamma[j,1:L[j]])))#
nx <- length(x[clust==k,1])
if (nx>0){
mm<-max(x[clust==k,1])#
curve(norm,from = 0, to = mm+0.5,col = k, lty=k,lwd=1, add = T)}}#
}#
dev.set()#
}#
#############################################################
}#
#
}#
}#
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################


#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
if (m==3){#
if (q==1){stop("The number of conditions is equal to 1. Choose m=1.")}#
print(paste("fitting the b_k parameterization"))#
print(paste("***************************************"))#
print(paste("*                K =", nc,"               *"))#
print(paste("***************************************"))#
run<-bkmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,t2,m2,prev.z,prev.clust,start.type=1,prev.alpha,prev.beta)#
icl<-run$icl#
bbs[nc,]<-c(run$bic,run$icl,run$ll)#
runs[[nc]]<-run#
for (nc in (Kmin+1):Kmax){#
print(paste("***************************************"))#
print(paste("*                K =", nc,"               *"))#
print(paste("***************************************"))#
z<-run$z#
ml<-length(run$psim)/(nc-1)#
alpha<-array(run$alpha[ml,,],dim=c(q,nc-1))#
beta<-array(run$beta[ml,,],dim=c(nc-1,tau))#
clust<-run$clust#
run <- tryCatch(
{pissa <- bkmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,tsplit,msplit,prev.z=z,prev.clust=clust,start.type=2,prev.alpha=alpha,prev.beta=beta)},
error = function(pissa){
	pissa <- bkmodel(reference=x,response=y,L,m=max.iter,K=nc,nr=-10*log(10),maxnr=maxnr,tsplit,msplit,prev.z=z,prev.clust=clust,start.type=1,prev.alpha=alpha,prev.beta=beta);return(pissa)
},
finally = pissa
)

bbs[nc,]<-c(run$bic,run$icl,run$ll)#
runs[[nc]]<-run#
plot(c(1,nc),c(min(bbs[1:nc,1],na.rm=T),max(bbs[1:nc,2],na.rm=T)),type="n",xlab="K",ylab="criterion")#
lines(bbs[1:nc,1])#
lines(bbs[1:nc,2],col=2)#
legend("topright",col=c(1,2),lty=1,c("BIC","ICL"))#
# graphics#
if (dn==T){#
if(run$icl<icl){#
icl<-run$icl#
print(paste("found new best model: K=",nc))#
dev.set()#
#############################################################
yy<-array(data = NA, dim = c(dim(y)[1],q))#
for (j in 1:q)yy[,j]<-rowSums(array(y[,index[j]:(index[j]+L[j]-1)],dim=c(dim(y)[1],L[j])))#
gamma<-run$gamma#
ml<-length(run$psim)/nc#
alpha<-array(run$alpha[ml,,],dim=c(q,nc))#
beta<-array(run$beta[ml,,],dim=c(nc,tau))#
op<-par(mfrow=c(q,2))#
clust<-run$clust#
for(j in 1:q){#
yyy<-log(yy[,j])#
plot(as.data.frame(cbind(x,exp(yyy))),col=run$clust,ylab=paste("sum(", L[j],"replicates)"),xlab="x",main=paste("condition: ", j))#
for (k in 1:nc){#
norm <- function(u) exp(alpha[j,k] + beta[k,1]*u + log(sum(exp(gamma[j,1:L[j]]))))#
nx <- length(x[clust==k,1])
if (nx>0){
mm<-max(x[clust==k,1])#
curve(norm,from = 0, to = mm+0.5,col = k, lty=k,lwd=1, add = T)}}#
#
plot(as.data.frame(cbind(x,yyy)),col=run$clust,ylab=paste("log(sum(", L[j],"replicates))"),xlab="x",main=paste("condition: ", j))#
for (k in 1:nc){#
norm <- function(u) alpha[j,k] + beta[k,1]*u + log(sum(exp(gamma[j,1:L[j]])))#
nx <- length(x[clust==k,1])
if (nx>0){
mm<-max(x[clust==k,1])#
curve(norm,from = 0, to = mm+0.5,col = k, lty=k,lwd=1, add = T)}}#
}#
dev.set()#
}#
#
#############################################################
}#
#
}#
#
}#
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################

colnames(bbs)<-c("BIC","ICL","Loglikelihood")
sel.mod.icl<-order(bbs[,2])[1]
sel.mod.bic<-order(bbs[,1])[1]


# final estimates for the selected model according to ICL criterion
r.icl <- runs[[sel.mod.icl]]
nc <- sel.mod.icl
tem <- length(r.icl$psim)/nc
est.sel.mod.icl <- vector("list",length=6)
if(nc == 1){
	est.sel.mod.icl[[1]] <- 1
}else{
	est.sel.mod.icl[[1]] <- r.icl$psim[tem,]
}
est.sel.mod.icl[[2]] <- r.icl$alpha[tem,,]
if (m==1){
est.sel.mod.icl[[3]] <- r.icl$beta[tem,,,]}
if (m==2){
est.sel.mod.icl[[3]] <- r.icl$beta[tem,,]}
if (m==3){
est.sel.mod.icl[[3]] <- r.icl$beta[tem,,]}
est.sel.mod.icl[[4]] <- r.icl$gamma
est.sel.mod.icl[[5]] <- r.icl$z
est.sel.mod.icl[[6]] <- r.icl$clust
names(est.sel.mod.icl)<-c("pi","alpha","beta","gamma","tau","clust")


# final estimates for the selected model according to BIC
r.bic <- runs[[sel.mod.bic]]
nc <- sel.mod.bic
tem <- length(r.bic$psim)/nc
est.sel.mod.bic <- vector("list",length=6)
if(nc == 1){
est.sel.mod.bic[[1]] <- 1
}else{
est.sel.mod.bic[[1]] <- r.bic$psim[tem,]
}
est.sel.mod.bic[[2]] <- r.bic$alpha[tem,,]
if (m==1){
est.sel.mod.bic[[3]] <- r.bic$beta[tem,,,]}
if (m==2){
est.sel.mod.bic[[3]] <- r.bic$beta[tem,,]}
if (m==3){
est.sel.mod.bic[[3]] <- r.bic$beta[tem,,]}
est.sel.mod.bic[[4]] <- r.bic$gamma
est.sel.mod.bic[[5]] <- r.bic$z
est.sel.mod.bic[[6]] <- r.bic$clust
names(est.sel.mod.bic)<-c("pi","alpha","beta","gamma","tau","clust")


results<-list(bbs,runs,sel.mod.icl,sel.mod.bic,est.sel.mod.icl,est.sel.mod.bic)
names(results)<-c("information.criteria","runs","sel.mod.icl","sel.mod.bic","est.sel.mod.icl","est.sel.mod.bic")
return(results)

}
