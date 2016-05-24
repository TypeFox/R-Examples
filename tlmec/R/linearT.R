################################################################################
############## Mixed-Effects Models With Censored Response #####################
################################################################################

## Student-t ##

## Linear Case ##

## EM Algorithm ##

################################################################################
#################### EM Algorithm / Censored Response ##########################
################################################################################

EMT<-function(cc,y,x,z,nj,nu,initial,criteria,diagnostic,error,iter.max){

print("The EM algorithm may take a long time to run ....")

GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)

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

  while(criterio > 0.001){

    count <- count + 1
    #print(count)
    soma1<- matrix(0,q1,q1)
    soma2<-0
    soma3<- matrix(0,p,p)
    soma4<- matrix(0,p,1)
    soma5<- matrix(0,p,p)
    
    tbi=matrix(0,m2,m)
    tui<-rep(0,m)
    tubi=matrix(0,m2,m)
    tubbi=matrix(0,m2,m2)
    tuybi=matrix(0,N,m2)
    tuyyi=matrix(0,N,N)
    tuyi=matrix(0,N,m)
    ver<-matrix(0,m,1)

    for (j in 1:m ){

      cc1=cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]

      y1=y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      x1=matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1=matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)

      gammai=x1%*%beta1
      Psi<-(sigmae*diag(nj[j])+(z1)%*%D1%*%t(z1))
      Psi<-(Psi+t(Psi))/2
      delta<- solve((sigmae*iD1+t(z1)%*%(z1)))
      dm<-t(y1-gammai)%*%solve(Psi)%*%(y1-gammai)
      cdm<-as.numeric((nu+nj[j])/(nu+dm))

      if(sum(cc1)==0){

        tuy<- (matrix(y1,nj[j],1))*cdm
        tuyy<- (y1%*%t(y1))*cdm
        tu<-cdm
        tub<- (delta%*%t(z1)%*%(tuy-gammai*cdm))
        tubb<- delta%*%(diag(q1)*sigmae+(t(z1)%*%(tuyy-tuy%*%t(gammai)-gammai%*%t(tuy)+cdm*gammai%*%t(gammai))%*%z1%*%delta))
        tuyb<- ((tuyy-tuy%*%t(gammai))%*%z1%*%delta)
        tb<- (delta%*%t(z1)%*%(matrix(y1,nj[j],1)-gammai))
        
        etuy2<- (y1%*%t(y1))*(cdm^2)
        etuy <- (matrix(y1,nj[j],1))*(cdm^2)
        etuy0<- (cdm^2)

        vari2<-etuy2-etuy%*%t(gammai)-gammai%*%t(etuy)+etuy0*(gammai%*%t(gammai))
        vari<-(tuy-(tu*gammai))%*%t(tuy-(tu*gammai))

        ver[j,]<- dmvt(as.vector(y1),as.vector(gammai),as.matrix(Psi),df=nu,log=FALSE)

      }


      if(sum(cc1)>0){

         if(sum(cc1)==nj[j]){

          muUi<-gammai
          SigmaUi<-Psi
          SigmaUiA<-SigmaUi*nu/(nu+2)
          auxupper<-y1-muUi
          auxU1<-pmvt(upper = c(auxupper), sigma = SigmaUiA, df = nu+2,algorithm = GB)[1]
          auxU2<-pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu,algorithm = GB) [1]
          MoMT = Mtmvt(muUi,SigmaUiA,nu+2,rep(-Inf,nj[j]),y1)
          MoMT1 = Mtmvt(muUi,SigmaUi,nu,rep(-Inf,nj[j]),y1)
          U0<-as.numeric(auxU1/auxU2)
          U1<-auxU1/auxU2*MoMT$Ey
          U2<-auxU1/auxU2*MoMT$Eyy

          tuy<-U1
          tuyy<-U2
          tu<-U0
          tub<-delta%*%t(z1)%*%(U1-gammai*U0)
          tubb<-delta%*%(diag(q1)*sigmae+(t(z1)%*%(U2-U1%*%t(gammai)-gammai%*%t(U1)+U0*gammai%*%t(gammai))%*%z1%*%delta))
          tuyb<-(U2-U1%*%t(gammai))%*%z1%*%delta
          tb<- (delta%*%t(z1)%*%(matrix(MoMT1$Ey,nj[j],1)-gammai))

          cp<-(((nu+nj[j])/nu)^2)*((gamma((nu+nj[j])/2)*gamma((nu+4)/2))/(gamma(nu/2)*gamma((nu+nj[j]+4)/2)))

          SigmaUiAey<-SigmaUi*nu/(nu+4)
          auxU1ey<-pmvt(upper = c(auxupper), sigma = SigmaUiAey, df = nu+4,algorithm = GB)[1]
          MoMTey = Mtmvt(muUi,SigmaUiAey,nu+4,rep(-Inf,nj[j]),y1)

          etuy2<- cp*(auxU1ey/auxU2)*MoMTey$Eyy
          etuy<- cp*(auxU1ey/auxU2)*MoMTey$Ey
          etuy0<- cp*(auxU1ey/auxU2)

          vari2<-etuy2-etuy%*%t(gammai)-gammai%*%t(etuy)+etuy0*(gammai%*%t(gammai))
          vari<-(tuy-(tu*gammai))%*%t(tuy-(tu*gammai))


          ver[j,]<-pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu,algorithm = GB) [1]

        }

        else {

          PsiA<-Psi*nu/(nu+2)
          nu1<-(nu+length(cc1[cc1==0]))

          muc<-x1[cc1==1,]%*%beta1+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)
          Sc <-Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]
          ScA<-nu/(nu+2)*Sc

          Qy1<-t(y1[cc1==0]-x1[cc1==0,]%*%beta1)%*%solve(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)
          Qy2<-t(y1[cc1==0]-x1[cc1==0,]%*%beta1)%*%solve(PsiA[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)

          auxcte<-as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
          auxcte1<-as.numeric((nu+2+Qy2)/(nu+2+length(cc1[cc1==0])))

          Sc22<-auxcte*Sc

          muUi<-muc
          SigmaUi<-Sc22

          SigmaUiA<- auxcte1*ScA
          auxupper<- y1[cc1==1]-muUi


          auxU1<-pmvt(upper = c(auxupper), sigma = SigmaUiA, df = nu1+2, algorithm = GB)[1]
          auxU2<-pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu1, algorithm = GB)[1]

          MoMT = Mtmvt(muUi,SigmaUiA,nu1+2,rep(-Inf,length(cc1[cc1==1])),y1[cc1==1])
          MoMT1 = Mtmvt(muUi,SigmaUi,nu1,rep(-Inf,length(cc1[cc1==1])),y1[cc1==1])
          
          U0<-as.numeric(auxU1/auxU2)/auxcte
          U1<-(U0)*(MoMT$Ey)
          U2<-(U0)*(MoMT$Eyy)

          Auxtuy<-(matrix(y1,nj[j],1))

          tuy<- Auxtuy*U0
          tuy[cc1==1]<- U1

          tuyy<-(Auxtuy%*%t(Auxtuy))

          AAx<-tuyy[cc1==0,cc1==0]*U0
          ABx<-Auxtuy[cc1==0]%*%t(U1)
          BAx<-t(ABx)
          BBx<-U2

          tuyy[cc1==0,cc1==0]<- AAx
          tuyy[cc1==0,cc1==1]<- ABx
          tuyy[cc1==1,cc1==0]<- BAx
          tuyy[cc1==1,cc1==1]<- BBx

          tu<-U0
          tub<-delta%*%t(z1)%*%(tuy-gammai*tu)
          tubb<-delta%*%(diag(q1)*sigmae+(t(z1)%*%(tuyy-tuy%*%t(gammai)-gammai%*%t(tuy)+tu*gammai%*%t(gammai))%*%z1%*%delta))
          tuyb<-(tuyy-tuy%*%t(gammai))%*%z1%*%delta
          
          eee<-matrix(y1,nj[j],1)
          eee[cc1==1]<-MoMT1$Ey
          tb<- (delta%*%t(z1)%*%(eee-gammai))

          SigmaUiAey<- as.numeric((nu+Qy1)/(nu+4+length(cc1[cc1==0])))*Sc
          SigmaUiAey<-(SigmaUiAey+t(SigmaUiAey))/2

          auxU1ey<-pmvt(upper = c(auxupper), sigma = SigmaUiAey, df = nu1+4, algorithm = GB)[1]

          MoMTey = Mtmvt(muUi,SigmaUiAey,nu1+4,rep(-Inf,length(cc1[cc1==1])),y1[cc1==1])

          dp<-((nu+nj[j])^2)*((gamma((nu+nj[j])/2)*gamma((nu+4+length(cc1[cc1==0]))/2))/
                  (gamma((nu+length(cc1[cc1==0]))/2)*gamma((nu+4+nj[j])/2)))

          U0ey<-as.numeric((auxU1ey/auxU2)*(dp/((nu+Qy1)^2)))
          U1ey<-(U0ey)*(MoMTey$Ey)
          U2ey<-(U0ey)*(MoMTey$Eyy)

          Auxtuyey<-(matrix(y1,nj[j],1))

          etuy2<-(Auxtuyey%*%t(Auxtuyey))

          AAxey<-etuy2[cc1==0,cc1==0]*U0ey
          ABxey<-Auxtuyey[cc1==0]%*%t(U1ey)
          BAxey<-t(ABxey)
          BBxey<-U2ey

          etuy2[cc1==0,cc1==0]<- AAxey
          etuy2[cc1==0,cc1==1]<- ABxey
          etuy2[cc1==1,cc1==0]<- BAxey
          etuy2[cc1==1,cc1==1]<- BBxey

          etuy<- Auxtuyey*U0ey
          etuy[cc1==1]<-U1ey

          etuy0<-U0ey

          vari2<-etuy2-etuy%*%t(gammai)-gammai%*%t(etuy)+etuy0*(gammai%*%t(gammai))
          vari<-(tuy-(tu*gammai))%*%t(tuy-(tu*gammai))

          SigmaUi<-(SigmaUi+t(SigmaUi))/2
          Psi[cc1==0,cc1==0]<-(Psi[cc1==0,cc1==0]+t(Psi[cc1==0,cc1==0]))/2

          ver[j,]<-dmvt(y1[cc1==0],gammai[cc1==0],as.matrix(Psi[cc1==0,cc1==0]),df=nu,log=FALSE)*(pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu1, algorithm = GB)[1])

        }

      }

      soma1<- soma1 + tubb
      soma2<- soma2 + (sum(diag(tuyy))-t(tuy)%*%gammai-t(gammai)%*%tuy-sum(diag(t(tuyb)%*%z1))-sum(diag(tuyb%*%t(z1)))
                          +t(gammai)%*%z1%*%tub+t(tub)%*%t(z1)%*%gammai+tu*(t(gammai)%*%gammai)+sum(diag(tubb%*%t(z1)%*%z1)))
      soma3<- soma3 + (tu*(t(x1)%*%x1))
      soma4<- soma4 + (t(x1)%*%(tuy-z1%*%tub))
      soma5<- soma5 + ((((nu+nj[j])/(nu+nj[j]+2))*t(x1)%*%solve(Psi)%*%x1)-(t(x1)%*%solve(Psi)%*%(vari2-vari)%*%solve(Psi)%*%x1))

      tubi[(((j-1)*q1)+1) : (j*q1), j]<-tub
      tubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]<-tubb
      tuybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]<-tuyb
      tuyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]<-tuyy
      tuyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]<-tuy
      tui[j]<-tu
      tbi[(((j-1)*q1)+1) : (j*q1), j]<-tb

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
     obj.out <- c(list(beta = beta1, sigmae= sigmae, dd = dd, varfixed = varbeta, bi = tbi, crit, iter = count, ubi = tubi, 
                  ubbi = tubbi, uybi = tuybi, uyi = tuyi, uyyi = tuyyi, ui = tui, diagnostic = diagnostic),crit)
  else
     obj.out <- c(list(beta = beta1, sigmae= sigmae, dd = dd, varfixed = varbeta, bi = tbi, loglik=loglik, AIC=AICc, BIC=BICc, 
                       AICcorr=AICcorr, iter = count, diagnostic = diagnostic),crit)

  class(obj.out) <- "t"

  return(obj.out)

}

