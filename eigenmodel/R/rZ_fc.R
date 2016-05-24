"rZ_fc" <-
function( EZ=XB(X,b)+ULU(UL), MH=TRUE ) {

## sample normal quantiles for probit model
## needs pp_zq, Ranks

  sd_zq<-1/sqrt(pp_zq)
  zq<-c(-Inf,rep(NA,max(Ranks,na.rm=TRUE)-1),Inf)
  for(ry in 1:(max(Ranks,na.rm=TRUE)-1)){
    ub<-suppressWarnings(min(Z[ Ranks==ry+1 ],na.rm=TRUE ) )
    lb<-suppressWarnings(max(Z[ Ranks==ry ],na.rm=TRUE ) )
    zq[ry+1]<-  qnorm( runif(1,pnorm(lb,0,sd_zq),pnorm(ub,0,sd_zq)),0,sd_zq  )
                                     }

## sample Z from fc  
## needs Ranks, uRanks to exist

  for(ry in sample(uRanks)){
    ir<- ( Ranks==ry & !is.na(Ranks) )
    lb<- zq[ry]
    ub<- zq[ry+1]
    z<-qnorm(runif(sum(ir),pnorm(lb,EZ[ir],1),pnorm(ub,EZ[ir],1)),EZ[ir],1)
    z[z== Inf]<-lb
    z[z==-Inf]<-ub
    Z[ir]<-z
                                  }
  ir<-is.na(Ranks)
  Z[ir]<-rnorm(sum(ir),EZ[ir],1)
  Z<-Z*upper.tri(Z)+ t(Z)*lower.tri(Z,diag=TRUE)
  diag(Z)<-NA


## MH proposal to help mixing
   if(MH) {
     del<-rnorm(1,0,1/sqrt(n))
     Zp<-Z+del ; zqp<-zq+del
     lhr<-sum(dnorm(Zp-EZ,0,1,log=TRUE),na.rm=TRUE)/2 -
          sum(dnorm(Z-EZ,0,1,log=TRUE),na.rm=TRUE)/2  +
          sum(dnorm(zqp[-c(1,length(uRanks)+1)],0,1/sqrt(pp_zq),log=TRUE)) -
          sum(dnorm(zq[-c(1,length(uRanks)+1)],0,1/sqrt(pp_zq),log=TRUE))
     if(log(runif(1))<lhr) { Z<-Zp ; zq<-zqp   }
          }

Z                        
                                            }

