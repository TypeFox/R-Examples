# x est le point etudie
#nbclust le nombre de classes sur xj
#mixmod le tableau mixmod pour xj
tik<-function(x=x,nbclust=nbclust,mixmod=mixmod){
   prob=rep(0,times=nbclust)
   for(i in 1:nbclust){
      prob[i]=mixmod[i,1]*dnorm(x,mean = mixmod[i,2],sd=sqrt(mixmod[i,3]))
      if(is.na(prob[i])){
         print("NA in tik")
         prob[i]=0
      }
   }
   sumprob=sum(prob)
   if(sumprob==0){
      prob[1]=1
   }else{
      prob=prob/sumprob
   }
   
   return(prob)
}