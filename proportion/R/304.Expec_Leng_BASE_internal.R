##### 1.WALD Expected Length for a given n and alpha level
gexplWD<-function(n,alp,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pW=0
  qW=0
  seW=0
  LW=0
  UW=0
  s=5000
  LEW=0

  ewiW=matrix(0,k,s)
  ewW=0
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #WALD METHOD
  for(i in 1:k)
  {
    pW[i]=x[i]/n
    qW[i]=1-(x[i]/n)
    seW[i]=sqrt(pW[i]*qW[i]/n)
    LW[i]=pW[i]-(cv*seW[i])
    UW[i]=pW[i]+(cv*seW[i])
    if(LW[i]<0) LW[i]=0
    if(UW[i]>1) UW[i]=1
    LEW[i]=UW[i]-LW[i]
  }
  #sumLEW=sum(LEW)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiW[i,j]=LEW[i]*dbinom(i-1, n,hp[j])
    }
    ewW[j]=sum(ewiW[,j])
  }
  ELW=data.frame(hp,ew=ewW,method="Wald")
  return(ELW)
}

##### 2.SCORE - Expected Length for a given n and alpha level
gexplSC<-function(n,alp,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pS=0
  qS=0
  seS=0
  LS=0
  US=0
  s=5000
  LES=0 								#LENGTH OF INTERVAL
  ewiS=matrix(0,k,s)						#Expected length quantity in sum
  ewS=0									#Expected Length
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  cv1=(cv^2)/(2*n)
  cv2=(cv/(2*n))^2

  #SCORE (WILSON) METHOD
  for(i in 1:k)
  {
    pS[i]=x[i]/n
    qS[i]=1-(x[i]/n)
    seS[i]=sqrt((pS[i]*qS[i]/n)+cv2)
    LS[i]=(n/(n+(cv)^2))*((pS[i]+cv1)-(cv*seS[i]))
    US[i]=(n/(n+(cv)^2))*((pS[i]+cv1)+(cv*seS[i]))
    if(LS[i]<0) LS[i]=0
    if(US[i]>1) US[i]=1
    LES[i]=US[i]-LS[i]
  }
  #sumLES=sum(LES)

  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiS[i,j]=LES[i]*dbinom(i-1, n,hp[j])
    }
    ewS[j]=sum(ewiS[,j])						#Expected Length
  }
  ELS=data.frame(hp,ew=ewS,method="Wilson")
  return(ELS)
}

##### 3. ARC SINE - Expected Length for a given n and alpha level
gexplAS<-function(n,alp,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pA=0
  qA=0
  seA=0
  LA=0
  UA=0
  s=5000
  LEA=0 								#LENGTH OF INTERVAL

  ewiA=matrix(0,k,s)						#Expected length quantity in sum
  ewA=0									#Expected Length
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #WALD METHOD
  #ARC-SINE METHOD
  for(i in 1:k)
  {
    pA[i]=x[i]/n
    qA[i]=1-pA[i]
    seA[i]=cv/sqrt(4*n)
    LA[i]=(sin(asin(sqrt(pA[i]))-seA[i]))^2
    UA[i]=(sin(asin(sqrt(pA[i]))+seA[i]))^2
    if(LA[i]<0) LA[i]=0
    if(UA[i]>1) UA[i]=1
    LEA[i]=UA[i]-LA[i]
  }
  #sumLEA=sum(LEA)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiA[i,j]=LEA[i]*dbinom(i-1, n,hp[j])
    }
    ewA[j]=sum(ewiA[,j])						#Expected Length
  }
  ELA=data.frame(hp,ew=ewA,method="ArcSine")
  return(ELA)
}


##### 4.LOGIT-WALD - Expected Length for a given n and alpha level
gexplLT<-function(n,alp,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pLT=0
  qLT=0
  seLT=0
  lgit=0
  LLT=0
  ULT=0
  s=5000
  LELT=0 								#LENGTH OF INTERVAL

  ewiLT=matrix(0,k,s)						#Expected length quantity in sum
  ewLT=0									#Expected Length
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #LOGIT-WALD METHOD

  pLT[1]=0
  qLT[1]=1
  LLT[1] = 0
  ULT[1] = 1-((alp/2)^(1/n))

  pLT[k]=1
  qLT[k]=0
  LLT[k]= (alp/2)^(1/n)
  ULT[k]=1

  for(j in 1:(k-2))
  {
    pLT[j+1]=x[j+1]/n
    qLT[j+1]=1-pLT[j+1]
    lgit[j+1]=log(pLT[j+1]/qLT[j+1])
    seLT[j+1]=sqrt(pLT[j+1]*qLT[j+1]*n)
    LLT[j+1]=1/(1+exp(-lgit[j+1]+(cv/seLT[j+1])))
    ULT[j+1]=1/(1+exp(-lgit[j+1]-(cv/seLT[j+1])))
    if(LLT[j+1]<0) LLT[j+1]=0
    if(ULT[j+1]>1) ULT[j+1]=1
  }
  for(i in 1:k)
  {
    LELT[i]=ULT[i]-LLT[i]
  }
  #sumLET=sum(LELT)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiLT[i,j]=LELT[i]*dbinom(i-1, n,hp[j])
    }
    ewLT[j]=sum(ewiLT[,j])						#Expected Length
  }
  ELLT=data.frame(hp,ew=ewLT,method="Logit-Wald")
  return(ELLT)
}

##### 5.t-WALD - Expected Length for a given n and alpha level
gexplTW<-function(n,alp,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pTW=0
  qTW=0
  seTW=0
  LTW=0
  UTW=0
  DOF=0
  cv=0
  s=5000
  LETW=0 								#LENGTH OF INTERVAL

  ewiTW=matrix(0,k,s)						#Expected length quantity in sum
  ewTW=0									#Expected Length
  #MODIFIED_t-WALD METHOD
  for(i in 1:k)
  {
    if(x[i]==0||x[i]==n)
    {
      pTW[i]=(x[i]+2)/(n+4)
      qTW[i]=1-pTW[i]
    }else
    {
      pTW[i]=x[i]/n
      qTW[i]=1-pTW[i]
    }
    f1=function(p,n) p*(1-p)/n
    f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
    DOF[i]=2*((f1(pTW[i],n))^2)/f2(pTW[i],n)
    cv[i]=qt(1-(alp/2), df=DOF[i])
    seTW[i]=cv[i]*sqrt(f1(pTW[i],n))
    LTW[i]=pTW[i]-(seTW[i])
    UTW[i]=pTW[i]+(seTW[i])
    if(LTW[i]<0) LTW[i]=0
    if(UTW[i]>1) UTW[i]=1
    LETW[i]=UTW[i]-LTW[i]
  }
  #sumLETW=sum(LETW)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiTW[i,j]=LETW[i]*dbinom(i-1, n,hp[j])
    }
    ewTW[j]=sum(ewiTW[,j])						#Expected Length
  }
  ELTW=data.frame(hp,ew=ewTW,method="Wald-T")
  # windows()
  # plot(ELTW,xlab="p",ylab="Expected Length",main="t-Wald",type="l")
  # abline(v=0.5, lty=2)
  return(ELTW)
}

#####6.LIKELIHOOD RATIO - Expected Length for a given n and alpha level
gexplLR<-function(n,alp,a,b)
{
  ####INPUT n
  y=0:n
  k=n+1
  ####INITIALIZATIONS
  mle=0
  cutoff=0
  LL=0
  UL=0
  s=5000
  LEL=0 								#LENGTH OF INTERVAL

  ewiL=matrix(0,k,s)						#Expected length quantity in sum
  ewL=0										#Simulation run to generate hypothetical p
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #LIKELIHOOD-RATIO METHOD
  for(i in 1:k)
  {
    likelhd = function(p) dbinom(y[i],n,p)
    loglik = function(p) dbinom(y[i],n,p,log=TRUE)
    mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
    cutoff[i]=loglik(mle[i])-(cv^2/2)
    loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
    LL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
    UL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
    LEL[i]=UL[i]-LL[i]
  }
  #sumLEL=sum(LEL)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiL[i,j]=LEL[i]*dbinom(i-1, n,hp[j])
    }
    ewL[j]=sum(ewiL[,j])						#Expected Length
  }
  ELL=data.frame(hp,ew=ewL,method="likelihood")
  return(ELL)
}

