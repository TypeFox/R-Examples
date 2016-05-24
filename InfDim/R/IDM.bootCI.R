IDM.bootCI <-
function(data,age,nboot){

#calculates covariance matrix from data and does the IDM calculations
P=cov(data)
out=IDM(P=P,age=age)

#Bootstrapping
Eig.vec=out$Eigenvectors
Eig.val=out$Eigenvalues
out.Traj=out$Trajectories
proc.Traj=Eig.val/sum(Eig.val)

ndata=length(data[,1])

boot.eigval=array(data=NA,dim=c(length(Eig.val),nboot))
boot.proc.traj=array(data=NA,dim=c(length(Eig.val),nboot))
boot.eigvec=array(data=NA,dim=c(nrow(Eig.vec),ncol(Eig.vec),nboot))
boot.traj=array(data=NA,dim=c(nrow(out.Traj),ncol(out.Traj),nboot))

for(i in 1:nboot){

ind.newdata=sample(1:ndata,replace=TRUE)
newdata=data[ind.newdata,]
Pnew=cov(newdata)

out.new=IDM(P=Pnew,age=age)

Eig.vec.new=out.new$Eigenvectors
Eig.val.new=out.new$Eigenvalues
out.Traj.new=out.new$Trajectories

for(nncol in 1:ncol(Eig.vec.new)){
cor.coef=cor(Eig.vec[,nncol],Eig.vec.new[,nncol])
if(cor.coef<0){
Eig.vec.new[,nncol]=-1*Eig.vec.new[,nncol]
out.Traj.new[,nncol]=-1*out.Traj.new[,nncol]
}
}

for(nn in 1:length(Eig.val.new)){
if(Eig.val.new[nn]<0){
value=Eig.val.new[nn]
Eig.val.new[nn]=0
Eig.val.new[1]=Eig.vec.new[1]-abs(value)
}}

boot.eigval[,i]=Eig.val.new
boot.proc.traj[,i]=Eig.val.new/sum(Eig.val.new)
boot.eigvec[,,i]=Eig.vec.new
boot.traj[,,i]=out.Traj.new
}

#CI for eigenvalues
eigval.ci.boot0.025=rep(NA,length(Eig.val))
eigval.ci.boot0.975=rep(NA,length(Eig.val))
proc.traj0.025=rep(NA,length(Eig.val))
proc.traj0.975=rep(NA,length(Eig.val))

for(i in 1:length(Eig.val)){
eigval.ci.boot0.025[i]=quantile(boot.eigval[i,],probs=c(0.025))
eigval.ci.boot0.975[i]=quantile(boot.eigval[i,],probs=c(0.975))
proc.traj0.025[i]=quantile(boot.proc.traj[i,],probs=c(0.025))
proc.traj0.975[i]=quantile(boot.proc.traj[i,],probs=c(0.975))
}

#CI for trajectories
traj.ci.boot0.025=matrix(NA,ncol=length(Eig.val),nrow=length(Eig.val))
traj.ci.boot0.975=matrix(NA,ncol=length(Eig.val),nrow=length(Eig.val))

for(i in 1:length(Eig.val)){
for(j in 1:length(Eig.val)){
traj.ci.boot0.025[i,j]=quantile(boot.traj[i,j,],probs=c(0.025))
traj.ci.boot0.975[i,j]=quantile(boot.traj[i,j,],probs=c(0.975))
}}

return(list(Eigenvalues=Eig.val,Eigenvalues2.5CI=eigval.ci.boot0.025, 
Eigenvalues97.5CI=eigval.ci.boot0.975,
Percent.trajectory=proc.Traj*100,
Percent.trajectory2.5CI=proc.traj0.025*100,
Percent.trajectory97.5CI=proc.traj0.975*100,
Trajectories=out.Traj,
Trajectories2.5CI=traj.ci.boot0.025,
Trajectories97.5CI=traj.ci.boot0.975))

}

