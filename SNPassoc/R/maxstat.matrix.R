maxstat.matrix<-function(x,...)
 {

########################################################################
#### auxiliary functions

      pmax<-function(t,rho,...){
      # Use abseps = 0.00001 to improve accuracy
        mean <- rep(0, 3)
        lower <- rep(-sqrt(t), 3)
        upper <- rep(sqrt(t), 3)
        corr <- diag(3)
        corr[1,2]<-corr[2,1]<-rho[1]
        corr[1,3]<-corr[3,1]<-rho[2]
        corr[2,3]<-corr[3,2]<-rho[3]
        prob <- pmvnorm(lower, upper, mean, corr=sqrt(corr), ...)
        prob
      }

      dGrec<-function(x) {
        aux<-prop.table(x,1)
        p<-aux[1,]
        q<-aux[2,]
        R<-apply(x,1,sum)[1]
        S<-apply(x,1,sum)[2]

        out<-rep(NA,6)
        out[1]<-out[2]<-((2* ((((q[1] + q[2]))* ((R - S)) + ((p[1] + p[2] + q[1] + q[2]))* R* log((2*((p[1] + p[2])))/(p[1] + p[2] + q[1] + q[2])))))/(p[1] + p[2] + q[1] + q[2]))
        out[3]<-((2* ((q[3]* ((R - S)) + ((p[3] + q[3]))* R*log((2* p[3])/(p[3] + q[3])))))/(p[3] + q[3]))
        out[4]<-out[5]<-((2* ((((p[1] + p[2]))* ((S - R)) + ((p[1] + p[2] + q[1] + q[2]))* S* log((2*((q[1] + q[2])))/(p[1] + p[2] + q[1] + q[2])))))/(p[1] + p[2] + q[1] + q[2]))
        out[6]<-((2* ((p[3]* ((S - R)) + ((p[3] + q[3]))* S*log((2* q[3])/(p[3] + q[3])))))/(p[3] + q[3]))

        out
      }

      dGdom<-function(x) {

        aux<-prop.table(x,1)
        p<-aux[1,]
        q<-aux[2,]
        R<-apply(x,1,sum)[1]
        S<-apply(x,1,sum)[2]

        out<-rep(NA,6)
        out[1]<- ((2*((q[1]* ((R - S)) + ((p[1] + q[1]))* R* log((2* p[1])/(p[1] + q[1])))))/(p[1] + q[1]))
        out[2]<-out[3]<-((2* ((((q[2] + q[3]))* ((R - S)) + ((p[2] + p[3] + q[2] + q[3])) *R* log((2*((p[2] + p[3])))/(p[2] + p[3] + q[2] + q[3])))))/(p[2] + p[3] + q[2] + q[3]))
        out[4]<- ((2*((p[1]* ((S - R)) + ((p[1] + q[1]))* S* log((2* q[1])/(p[1] + q[1])))))/(p[1] + q[1]))
        out[5]<-out[6]<-((2* ((((p[2] + p[3]))* ((S - R)) + ((p[2] + p[3] + q[2] + q[3])) *S* log((2*((q[2] + q[3])))/(p[2] + p[3] + q[2] + q[3])))))/(p[2] + p[3] + q[2] + q[3]))

        out
      }

      dGadd<-function(x) {

        aux<-prop.table(x,1)
        p<-aux[1,]
        q<-aux[2,]
        R<-apply(x,1,sum)[1]
        S<-apply(x,1,sum)[2]

        out<-rep(NA,6)

        out[1]<-((-(((((p[1] - p[3] - q[1] + q[3]))* R* S *((R + S))^2* (((((- q[1]) + p[1]* (((-1) + 2* q[1] - 2* q[3])) + q[3] + p[3]* (((-3) - 2 *
                  q[1] + 2* q[3]))))* R + 2* ((q[1]^2 + (((-1) + q[3]))* q[3] - q[1]* ((1 + 2* q[3]))))* S))))/((
                  p[1]^2* R^2 + (((-1) + p[3])) * p[3]* R^2 - ((p[3] + q[1] + 2* p[3]* q[1] +
                  q[3] - 2* p[3]* q[3]))* R * S + ((q[1]^2 + (((-1) + q[3])) *q[3] - q[1]* ((1 + 2 *
                  q[3])))) * S^2 - p[1]* R* ((R + 2* p[3]* R + S - 2 *q[1] * S + 2* q[3]* S))))^2)))
        out[2]<-0
        out[3]<- (((((p[1] - p[3] - q[1] + q[3]))* R* S* ((R + S))^2* ((((
                  q[1] + p[1]* (((-3) + 2* q[1] - 2* q[3])) - q[3] + p[3]* (((-1) - 2* q[1]
                  + 2* q[3]))))* R + 2* ((q[1]^2 + (((-1) + q[3]))* q[3] - q[1]* ((1 + 2* q[3]))))* S))))/((
                  p[1]^2* R^2 + (((-1) + p[3]))* p[3]* R^2 - ((p[3] + q[1] + 2 *p[3]* q[1] +
                  q[3] - 2 *p[3] * q[3]))* R* S + ((q[1]^2 + (((-1) + q[3])) * q[3] - q[1]* ((1 + 2* q[3]))))* S^2 -
                  p[1]* R* ((R + 2* p[3]* R + S - 2* q[1]* S + 2 * q[3] *S))))^2)
        out[4]<-(((( p[1] - p[3] - q[1] + q[3]))* R* S* ((R + S))^2*((
                  2* ((p[1]^2 + (((-1) + p[3]))* p[3] - p[1]* ((1 + 2* p[3]))))* R +
                  (((-p[1]) + p[3] - q[1] + 2* p[1]* q[1] - 2* p[3]* q[1] - 3 *q[3] - 2* p[1]* q[3] +
                  2* p[3]* q[3]))* S))))/((p[1]^2* R^2 + (((-1) + p[3]))* p[3]* R^2 - ((p[3] + q[1] + 2* p[3] *q[1] + q[3] - 2 *
                  p[3]* q[3]))* R* S + ((q[1]^2 + (((-1) + q[3]))* q[3] - q[1]* ((1 + 2* q[3])))) *S^2 - p[1]* R* ((R + 2* p[3]* R +
                  S - 2* q[1]* S + 2* q[3] *S))))^2
        out[5]<- 0
        out[6]<- (-(((((p[1] - p[3] - q[1] + q[3]))* R* S* ((R + S))^2 *
                 ((2* p[1]^2* R + 2* (((-1) + p[3]))* p[3]* R - ((p[3] + 3* q[1] + 2* p[3]* q[1] + q[3]
                 - 2* p[3]* q[3]))* S + p[1]* (((-2)* R - 4* p[3]* R + S + 2* q[1]* S - 2* q[3]* S)))
                 )))/((p[1]^2* R^2 + (((-1) + p[3]))* p[3]* R^2 - ((p[3] + q[1] + 2* p[3]* q[1] +
                 q[3] - 2* p[3]* q[3]))* R *S + ((q[1]^2 + (((-1) + q[3]))* q[3] - q[1]* ((1 + 2*
                 q[3])))) *S^2 - p[1]* R* ((R + 2* p[3]* R + S - 2* q[1]* S + 2* q[3] * S))))^2))

        out
      }

      sigmaMatrix<-function(x)
       {
        aux<-prop.table(x,1)
        p<-aux[1,]
        q<-aux[2,]
        ansp<-(diag(p)-(p%*%t(p)))
        ansq<-(diag(q)-(q%*%t(q)))

        ans<-matrix(0,nrow=6,ncol=6)
        ans[1:3,1:3]<-ansp
        ans[4:6,4:6]<-ansq

        ans
       }

########################################################################

# x is a 2x3 table where rows are case/control
  if(nrow(x)==1 | ncol(x)==1) 
    {
     ans<-list(stat=NA,pval=NA)
     tests<-rep(NA,5)
     names(tests)<-c("dominant","recessive","log-additive","MAX-statistic", "Pr(>z)")
    }
  else 
   {
    x[x==0]<-0.5
  
    if(ncol(x)==2)
     {
      estad<-G(x)
      dom<-estad
      rec<-estad
      add<-estad
      stat<-estad
      pval<-1-pchisq(stat,1) 
      ans<-list(stat=stat,pval=pval)
      warning("this SNP only has 2 genotypes")
     }
    else
     {
     # dominant
      xd<-matrix(NA,ncol=2,nrow=2)
      xd[,1]<-x[,1]
      xd[,2]<-x[,2]+x[,3]
      dom<-G(xd)

     # recessive
      xr<-matrix(NA,ncol=2,nrow=2)
      xr[,1]<-x[,1]+x[,2]
      xr[,2]<-x[,3]
      rec<-G(xr)

     # additive 
      add<-CA(x)

     # derivatives
      dd<-dGdom(x)
      rr<-dGrec(x)
      ad<-dGadd(x)

     # To control G-statistic equal to 0
     if (dom==0)
      dd<-dd+runif(length(dd))/1000 
     if (rec==0)
      rr<-rr+runif(length(rr))/1000
     if (add==0)
      ad<-ad+runif(length(ad))/1000


     # sigma
      ss<-sigmaMatrix(x)

     # correlation matrix
      rho<-rep(NA,3)
      rho[1]<-(t(dd)%*%ss%*%rr/(((t(dd)%*%ss%*%dd)^.5)*(t(rr)%*%ss%*%rr)^.5))^2
      rho[2]<-(t(dd)%*%ss%*%ad/(((t(dd)%*%ss%*%dd)^.5)*(t(ad)%*%ss%*%ad)^.5))^2
      rho[3]<-(t(rr)%*%ss%*%ad/(((t(rr)%*%ss%*%rr)^.5)*(t(ad)%*%ss%*%ad)^.5))^2

     # p-value
      stat<-max(max(dom,rec,add))
      pval<-1-pmax(stat,rho=rho,...)
      ans<-list(stat=stat,pval=pval)
     }
     tests<-c(dom,rec,add,stat, pval)
     names(tests)<-c("dominant","recessive","log-additive","MAX-statistic", "Pr(>z)")
    }
  class(tests)<-"maxstat"
  attr(tests,"info")<-ans
  tests
}




G<-function(tt)
{
 df<-(dim(tt)[1]-1)*(dim(tt)[2]-1)
 tt.r<-apply(tt,1,sum)
 R<-tt.r[1]
 S<-tt.r[2]
 N<-R+S
 n<-apply(tt,2,sum)
 a1<-(tt[1,]*N)/(R*n)
 a2<-(tt[2,]*N)/(S*n)
 ans <- 2*sum(tt*log(rbind(a1,a2)))
 ans
}



CA<-function(x)
 {
  aux<-prop.table(x,1)
  p<-aux[1,]
  q<-aux[2,]
  R<-apply(x,1,sum)[1] 
  S<-apply(x,1,sum)[2]
  N<-R+S

  num<-(p[1] - p[3] - q[1] + q[3])^2 * R* S* (R + S)
  den<- -((p[1] - p[3])* R + (q[1] - q[3])* S)^2 + (R + S)* ((p[1] + p[3])*R
        + (q[1] + q[3])* S)

  ans<-num/den
  ans
 }

