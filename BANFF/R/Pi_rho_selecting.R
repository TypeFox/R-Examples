#####Using likelyhood for choosing pi and rho
Pi_rho_selecting=function(znew,rstat,piall,rhoall)
{
  mylog<-0
  like<-0
  for(i in 1:(length(piall)*sum(seq(length(rhoall))-1))){
    mu0=mean(rstat[which(znew[i,]==0)])
    mu1=mean(rstat[which(znew[i,]==1)])
    var0=var(rstat[which(znew[i,]==0)])
    var1=var(rstat[which(znew[i,]==1)])
    for(num1 in 1:length(rstat)){
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
  if (length(choice)!=0){choice=sample(choice,1)}
  return(choice)
}