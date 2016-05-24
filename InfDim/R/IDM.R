IDM <-
function(P,age){

adjusted_age=adjust(age)
psiiB=psii(adjusted_age)

inv_psii=solve(psiiB)
inv_tpsii=solve(t(psiiB))

CP=inv_psii%*%P%*%inv_tpsii

eig_first=eigen(CP)

Eigval=eig_first$values
Eigvec=eig_first$vectors

traj=matrix(NA,nrow=length(age),ncol=length(age))

for(i in 1:length(age)){
for(j in 1:length(age)){
traj[j,i]=sum(psiiB[j,]*Eigvec[,i])
}}

for(i in 1:length(age)){
traj[,i]=traj[,i]/(sqrt(sum(traj[,i]^2)))
}

perc.var=Eigval/sum(Eigval)

return(list(Eigenvalues=Eigval,Eigenvectors=Eigvec,Percent.trajectory=perc.var*100,Trajectories=traj))
}

