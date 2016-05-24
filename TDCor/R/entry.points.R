entry.points <-
function(rd_sub_norm,net_f,thr_bool_EP,MinTarNumber,MinProp,MaxEPNumber)

{

net_f=sign(net_f)



############# Predicting steady states



S0=to.bool(rd_sub_norm,thr=thr_bool_EP)[,1] # Observed state at t=0 (First time point)

S1=to.bool(rd_sub_norm,thr=thr_bool_EP)[,2] # Observed state at t=6h (Second time point)



SS=rep(0,dim(rd_sub_norm)[1])  # Vector of predicted steady states: 1 if positive regulator on and negative off else 0

iin=which(rowSums(abs(rd_sub_norm))>0)

for (i in iin)

{

posreg=which(net_f[i,]>0)

negreg=which(net_f[i,]<0)



if (length(posreg)>0 & length(negreg)>0)

{

if (sum(S0[posreg])>0 & sum(S0[negreg])==0)

{

SS[i]=1

}

}



if (length(posreg)>0 & length(negreg)==0)

{

if (sum(S0[posreg])>0)

{

SS[i]=1

}

}

}



# Comparison of the predicted SS and the observed S0: 0 if steady state, 1 if not



SSC=abs(SS-S0)

SSC[which(rowSums(abs(net_f))==0)]=0



############# Looking for the Master Regulator-Signal Transducer (MRST)



regpot=colSums(abs(net_f)*SSC)



# Filter on the minimum number of targets not at steady state



regpot[regpot<MinTarNumber]=0



# Filter on the high expression of the regulator



regpot=regpot*S0



# Filter on the proportion of targets that are not at steady state at t=0



rin=which(regpot!=0)

for (r in rin)

{

targets=which(net_f[,r]!=0)

if (sum(SSC[targets])/sum(abs(net_f)[,r])<MinProp)

{regpot[r]=0}

}



# Filter on the maximum number of entry points



if (length(which(regpot!=0))>MaxEPNumber)

{

min_score=rev(sort(regpot))[min(MaxEPNumber,length(regpot))]

regpot[regpot<min_score]=0

}

regpot[regpot>0]=1



# Changing state between first and second time point (used in df_TDCOR for pruning)



CS=S0-S1



output=list(EP=regpot,SSC=SSC,CS=CS)

}
