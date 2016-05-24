######Parameter Selection
HyperPara.Select=function(net,pvalue,piall,rhoall,n=30)
{
  rstat=Transfer(pvalue)
  wholeindex=Kmeans(rstat)
  pirhopair=Pi_rho_gen(piall,rhoall)
  znew=Generating_Zi(net,pirhopair,wholeindex,n)
  mylog<-0
  like<-0
  for(i in 1:(length(piall)*sum(seq(length(rhoall))-1))){
    for(num1 in 1:length(rstat)){
      
      mu0=mean(rstat[which(znew[i,]==0)])
      mu1=mean(rstat[which(znew[i,]==1)])
      var0=var(rstat[which(znew[i,]==0)])
      var1=var(rstat[which(znew[i,]==1)])
      if(znew[i,num1]==0){
        mylog[num1]<--(rstat[num1]-mu0)^2/(2*var0)-log(sqrt(var0))
      }else{
        mylog[num1]<--(rstat[num1]-mu1)^2/(2*var1)-log(sqrt(var1))
      }
      
    }
    like=c(like,sum(mylog))
    likenew=like[-1]
  }
  maxlike=max(likenew[!is.na(likenew)])
  choice=which(likenew==maxlike)
  mychoice=list()
  mychoice$pi0=pirhopair$pi0[choice]
  mychoice$rho0=pirhopair$rho0[choice]
  mychoice$rho1=pirhopair$rho1[choice]
  return(mychoice)
}
