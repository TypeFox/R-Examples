################################################################################
############## Mixed-Effects Models With Censored Response #####################
################################################################################

## Normal ##

## Linear Case ##

## EM Algorithm ##

################################################################################
#################### EM Algorithm / Censored Response ##########################
################################################################################

EMN<-function(cc,y,x,z,nj,initial,criteria,diagnostic,error,iter.max){

  m<-length(nj)[1]
  N<-sum(nj)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  m1<-m*p
  m2<-m*q1
  
  
  #valores iniciais
  beta1<-matrix(initial$beta,p,1)
  sigmae<- initial$sigma2
  D1<-initial$Delta
  iD1<-solve(D1)

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)])

  criterio<-1
  count<-0

  while(criterio > error){

    count <- count + 1
    #print(count)
    soma1<- matrix(0,q1,q1)
    soma2<-0
    soma3<- matrix(0,p,p)
    soma4<- matrix(0,p,1)
    soma5<- matrix(0,p,p)

    ubi=matrix(0,m2,m)
    ubbi=matrix(0,m2,m2)
    uybi=matrix(0,N,m2)
    uyyi=matrix(0,N,N)
    uyi=matrix(0,N,m)

    ver<-matrix(0,m,1)

    for (j in 1:m ){
 
      cc1=cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      y1=y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      x1=matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1=matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
      gammai=x1%*%beta1

      if(sum(cc1)==0){

        delta<- solve(iD1+(t(z1)%*%((z1*(1/sigmae)))))
        uy<- matrix(y1,nj[j],1)
        uyy<- y1%*%t(y1)
        ub<- delta%*%(t(z1)*(1/sigmae))%*%(uy-gammai)
        ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%z1%*%delta)
        uyb<- (uyy-uy%*%t(gammai))%*%(z1*(1/sigmae))%*%delta
        
        Psi<-(sigmae*diag(nj[j])+(z1)%*%D1%*%t(z1))
        Psi<-(Psi+t(Psi))/2
        ver[j,]<- dmvnorm(as.vector(y1),gammai,Psi)

      }

      if(sum(cc1)>=1){

        Psi<-(sigmae*diag(nj[j])+(z1)%*%D1%*%t(z1))
        Psi<-(Psi+t(Psi))/2

        if(sum(cc1)==nj[j]){
        
          muc=as.vector(x1%*%beta1)
          Sc<-as.matrix(Psi)
          delta<- solve(iD1+(t(z1)%*%((z1*(1/sigmae)))))
          aux<- MomemNT(muc,Sc,y1)
          uy<-aux$Ey
          uyy<- aux$Eyy
          ub<- delta%*%(t(z1)*(1/sigmae))%*%(uy-gammai)
          ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%z1%*%delta)
          uyb<- (uyy-uy%*%t(gammai))%*%(z1*(1/sigmae))%*%delta

          ver[j,]<-pmvnorm(upper=y1,mean=muc,sigma=Sc)

        }

        else {
        
          muc=as.vector(x1[cc1==1,]%*%beta1+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1))
          Sc <-Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]
          delta<- solve(iD1+(t(z1)%*%((z1*(1/sigmae)))))
          Sc<-(Sc+t(Sc))/2
          dimnames(Sc) <- NULL
          aux <-MomemNT(muc,Sc,y1[cc1==1])
          uy <-matrix(y1,nj[j],1)
          uy[cc1==1]<-aux$Ey

          uyy<-matrix(0,nj[j],nj[j])
          uyy[cc1==1,cc1==1]<-aux$Vary
          uyy<- uyy+uy%*%t(uy)
          ub<- delta%*%(t(z1)*(1/sigmae))%*%(uy-gammai)
          ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%z1%*%delta)
          uyb<- (uyy-uy%*%t(gammai))%*%(z1*(1/sigmae))%*%delta

          ver[j,]<-dmvnorm(y1[cc1==0],gammai[cc1==0],as.matrix(Psi[cc1==0,cc1==0]))*pmvnorm(upper=y1[cc1==1],mean=muc,sigma=Sc)

        }

      }

      soma1<- soma1 + ubb
      soma2<- soma2 + (sum(diag(uyy))-t(uy)%*%gammai-t(gammai)%*%uy-sum(diag(t(uyb)%*%z1))-sum(diag(uyb%*%t(z1)))
                          +t(gammai)%*%z1%*%ub+t(ub)%*%t(z1)%*%gammai+t(gammai)%*%gammai+sum(diag(ubb%*%t(z1)%*%z1)))
      soma3<- soma3 + (t(x1)%*%x1)
      soma4<- soma4 + (t(x1)%*%(uy-z1%*%ub))
      soma5<- soma5 + (t(x1)%*%solve(Psi)%*%x1-t(x1)%*%solve(Psi)%*%(uyy-uy%*%t(uy))%*%solve(Psi)%*%x1)

      ubi[(((j-1)*q1)+1) : (j*q1), j]<-ub
      ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]<-ubb
      uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]<-uyb
      uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]<-uyy
      uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]<-uy

    }

    beta1<- solve(soma3)%*%soma4
    sigmae<- (1/N)*(soma2)
    sigmae<-as.numeric(sigmae)
    D1<- (1/m)*(soma1)
    iD1<- solve(D1)
    teta1 <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)])
    logver <- sum(log(ver))
    varbeta<-solve(soma5)

    if (count>1){
      criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
    }

    if (count==iter.max){
      criterio <- -1
    }

    teta<-teta1
    logver1<-logver

  }

  dd<-D1[upper.tri(D1, diag = T)]

  npar<-length(c(teta1))

  ni<-sum(nj)
  crit <- NULL
  if(criteria){
    loglik<-logver1
    AICc<- -2*loglik +2*npar
    AICcorr<- AICc + ((2*npar*(npar+1))/(ni-npar-1))
    BICc <- -2*loglik +log(ni)*npar
    crit <- list(loglik=loglik, AIC=AICc, BIC=BICc, AICcorr=AICcorr)
  }
  if(diagnostic)
     obj.out <- c(list(beta = beta1, sigmae= sigmae, dd = dd, varfixed = varbeta, bi = ubi,crit, iter = count, ubi = ubi, ubbi = ubbi, 
                  uybi = uybi, uyi = uyi, uyyi = uyyi, diagnostic = diagnostic),crit)
  else
     obj.out <- c(list(beta = beta1, sigmae= sigmae, dd = dd, varfixed = varbeta, bi = ubi,loglik=loglik, AIC=AICc, BIC=BICc, 
                       AICcorr=AICcorr, iter = count, diagnostic = diagnostic),crit)

  class(obj.out) <- "Normal"

  return(obj.out)
  
}

