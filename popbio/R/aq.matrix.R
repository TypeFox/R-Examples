aq.matrix<-function(trans, recruits, summary=TRUE, seed.survival=0.126, seed.bank.size=10000, seeds.per.fruit=120, ...)
{
   x<-trans
   seeds.from.plants<-sum(x$fruits)*seeds.per.fruit
   ## assume seeds in seed bank and new seeds have equal chance for successful germination
   recruitment.rate <- recruits/(seed.bank.size + seeds.from.plants)

   ## add fertilities
   x<-cbind(x, recruit=x$fruits / sum(x$fruits) * seeds.from.plants * recruitment.rate)
   x<-cbind(x, seed = x$fruits * seeds.per.fruit * seed.survival)
   if(summary)
   {
      ## STAGE vector
      n<- summary(x$stage)
      n["seed"]<- seed.bank.size
      ## matrix
    
      A <- projection.matrix(x, add=c(1,1,seed.survival, 2,1, recruitment.rate), ...)
      lam <- max(Re(eigen(A)$values))      

      n1 <- A %*% n
      ## format same as n
      n1<-as.vector(n1)
      names(n1)<-names(n)
    
      z<-list(
            recruits=recruits,
            seed.survival=seed.survival,
            seed.bank=seed.bank.size,
            seeds.from.plants=seeds.from.plants,
            recruitment.rate=recruitment.rate,
            A=A,
            lambda=lam,
            n=n,
            n1=round(n1,0)
              )
      z
   }
   else {x}
}

