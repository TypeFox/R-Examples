FPTime <-
function(state,nchains,tmat,io,n)
{
dm<-MultDTMC(nchains,tmat,io,n)
fp1<-c()
for(i in 1:nchains)
{
 dm2<-dm[[i]]
fp<-ifelse(dm2[state,]==1,1,0)
fp3<-which(fp>0)
#fp1[i]<-ifelse(fp3[1]==1,fp3[2],fp3[1]) use if we want first passage time after the first state
fp1[i]<-fp3[1]
}
return(fp1)
}

