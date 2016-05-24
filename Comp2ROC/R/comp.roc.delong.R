comp.roc.delong <-
function(sim1.ind,sim1.sta,sim2.ind,sim2.sta,related=TRUE) {
  
imax=length(sim1.sta[sim1.sta==0])
jmax=length(sim1.sta[sim1.sta==1])
imax2=length(sim2.sta[sim2.sta==0])
jmax2=length(sim2.sta[sim2.sta==1])

total_cases=imax+jmax
total_cases2=imax2+jmax2

dados=rbind(cbind(sim1.ind,sim1.sta),cbind(sim2.ind,sim2.sta))

mod=2

if(imax>imax2) maxImax = imax else maxImax = imax2
if(jmax>jmax2) maxJmax = jmax else maxJmax = jmax2

# norm defines the negative cases
norm=array(data=NA,c(maxImax,mod))
# abnorm defines the positive cases
abnorm=array(data=NA,c(maxJmax,mod))
# data are atributed to the arrays norm and abnorm
for (i in 1:imax) {
  norm[i,1]=dados[i,1]
}
for (j in 1:jmax) {
  abnorm[j,1]=dados[j+imax,1]
}
for (i in 1:imax2) {
  norm[i,2]=dados[i+total_cases,1]
}
for (j in 1:jmax2) {
  abnorm[j,2]=dados[j+total_cases+imax2,1]
}

# T1 and T2 defines the Wilcoxon Mann Whitney matrix for each modality
T1=array(data=NA,c(jmax,imax))
for (i in 1:imax){
  for (j in 1:jmax){
    dif=abnorm[j,1]-norm[i,1]
    if(dif>0) T1[j,i]=1 else (if(dif==0) T1[j,i]=0.5 else T1[j,i]=0)
  }
}

T2=array(data=NA,c(jmax2,imax2))
for (i in 1:imax2){
  for (j in 1:jmax2){
    dif=abnorm[j,2]-norm[i,2]
    if(dif>0) T2[j,i]=1 else (if(dif==0) T2[j,i]=0.5 else T2[j,i]=0)
  }
}

# U and V are the arrays obtained by summing the columns and the rows of T1
U=array(data=NA,c(maxJmax,mod))
V=array(data=NA,c(maxImax,mod))
for (j in 1:jmax){
  U[j,1]=sum(T1[j,])/imax
}
for (i in 1:imax){
  V[i,1]=sum(T1[,i])/jmax
}
for (j in 1:jmax2){
  U[j,2]=sum(T2[j,])/imax2
}
for (i in 1:imax2){
  V[i,2]=sum(T2[,i])/jmax2
}
# Areas under the ROC curves for each modality
AUC=array(data=NA,c(mod))
AUC[1]=1/(imax*jmax)*sum(T1[,])
AUC[2]=1/(imax2*jmax2)*sum(T2[,])

# Standard deviations for eachj modality
S=array(data=NA,c(mod))

S[1]=(1/imax*sum((V[,1]-AUC[1])^2, na.rm=TRUE)/(imax-1)+1/jmax*sum((U[,1]-AUC[1])^2, na.rm=TRUE)/(jmax-1))^(0.5)
S[2]=(1/imax2*sum((V[,2]-AUC[2])^2, na.rm=TRUE)/(imax2-1)+1/jmax2*sum((U[,2]-AUC[2])^2, na.rm=TRUE)/(jmax2-1))^(0.5)

# Vari?nces for each modality, positive and negative cases and global
S2V=array(data=NA,c(mod,mod))
S2U=array(data=NA,c(mod,mod))
S2=array(data=NA,c(mod,mod))

for (z in 1:mod){
  S2V[1,z]= sum((V[,1]-AUC[1])*(V[,z]-AUC[z]), na.rm=TRUE)/(imax-1)
  S2U[1,z]= sum((U[,1]-AUC[1])*(U[,z]-AUC[z]), na.rm=TRUE)/(jmax-1)
  S2[1,z]=1/imax*S2V[1,z]+1/jmax*S2U[1,z]
}
for (z in 1:mod){
  S2V[2,z]= sum((V[,2]-AUC[2])*(V[,z]-AUC[z]), na.rm=TRUE)/(imax2-1)
  S2U[2,z]= sum((U[,2]-AUC[2])*(U[,z]-AUC[z]), na.rm=TRUE)/(jmax2-1)
  S2[2,z]=1/imax2*S2V[2,z]+1/jmax2*S2U[2,z]
}

# Global correlations
R=array(data=NA,c(mod,mod))
for (k in 1:mod){
  for (z in 1:mod){
    if (related == TRUE) {
       if (z==k) {
        if (S2[k,k]*S2[z,z]==0)
          R[k,z]=1
        else
          R[k,z]=S2[k,z]/((S2[k,k]*S2[z,z])^(0.5))
          }
        else
        {
          if (S2[k,k]*S2[z,z]==0)
            R[k,z]=0
          else
            R[k,z]=S2[k,z]/((S2[k,k]*S2[z,z])^(0.5))
        }
    } else {
      R[k,z] = 0
    }
  }
}
Z=(AUC[1]-AUC[2])/sqrt(S2[1,1]+S2[2,2]-2*R[1,2]*sqrt(S2[1,1]*S2[2,2]))
if (Z>0) p.value=2*(1-pnorm(Z,0,1)) else p.value=2*pnorm(Z,0,1)

answer=list(Z=Z,pvalue=p.value,AUC=AUC,SE=S,S=S2,R=R)
return(answer)
}
