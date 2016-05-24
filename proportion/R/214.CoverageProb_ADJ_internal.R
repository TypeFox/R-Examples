##### 1.ADJUSTED WALD-Coverage Probability
gcovpAWD<-function(n,alp,h,a,b,t1,t2)
{
  ####INPUT n
  x=0:n
  k=n+1
  y=x+h
  n1=n+(2*h)
  ####INITIALIZATIONS
  pAW=0
  qAW=0
  seAW=0
  LAW=0
  UAW=0
  s=5000								#Simulation run to generate hypothetical p
  cpAW=matrix(0,k,s)
  ctAW=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppAW=0								#Coverage probabilty
  ctr=0
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #WALD METHOD
  for(i in 1:k)
  {
    pAW[i]=y[i]/n1
    qAW[i]=1-pAW[i]
    seAW[i]=sqrt(pAW[i]*qAW[i]/n1)
    LAW[i]=pAW[i]-(cv*seAW[i])
    UAW[i]=pAW[i]+(cv*seAW[i])
    if(LAW[i]<0) LAW[i]=0
    if(UAW[i]>1) UAW[i]=1
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LAW[i] && hp[j] < UAW[i])
      {
        cpAW[i,j]=dbinom(i-1, n,hp[j])
        ctAW[i,j]=1
      }
    }
    cppAW[j]=sum(cpAW[,j])
    if(t1<cppAW[j]&&cppAW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPAW=data.frame(hp,cp=cppAW,method="Adj-Wald")
  return(CPAW)
}

##### 2.ADJUSTED SCORE - Coverage Probability
gcovpASC<-function(n,alp,h,a,b,t1,t2)
{
  ####INPUT n
  x=0:n
  k=n+1
  y=x+h
  n1=n+(2*h)
  ####INITIALIZATIONS
  pAS=0
  qAS=0
  seAS=0
  LAS=0
  UAS=0
  s=1000								#Simulation run to generate hypothetical p
  cpAS=matrix(0,k,s)
  ctAS=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppAS=0							#Coverage probabilty
  ctr=0

  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  cv1=(cv^2)/(2*n1)
  cv2=(cv/(2*n1))^2

  #SCORE (WILSON) METHOD
  for(i in 1:k)
  {
    pAS[i]=y[i]/n1
    qAS[i]=1-pAS[i]
    seAS[i]=sqrt((pAS[i]*qAS[i]/n1)+cv2)
    LAS[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)-(cv*seAS[i]))
    UAS[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)+(cv*seAS[i]))
    if(LAS[i]<0) LAS[i]=0
    if(UAS[i]>1) UAS[i]=1
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LAS[i] && hp[j] < UAS[i])
      {
        cpAS[i,j]=dbinom(i-1, n,hp[j])
        ctAS[i,j]=1
      }
    }
    cppAS[j]=sum(cpAS[,j])						#Coverage Probability
    if(t1<cppAS[j]&&cppAS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPAS=data.frame(hp,cp=cppAS,method="Adj-Wilson")
  return(CPAS)
}


##### 3.ADJUSTED ARC SINE - Coverage Probability
gcovpAAS<-function(n,alp,h,a,b,t1,t2)
{
  ####INPUT n
  x=0:n
  k=n+1
  y=x+h
  n1=n+(2*h)
  ####INITIALIZATIONS
  pAA=0
  qAA=0
  seAA=0
  LAA=0
  UAA=0
  s=5000								#Simulation run to generate hypothetical p
  cpAA=matrix(0,k,s)
  ctAA=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppAA=0								#Coverage probabilty
  ctr=0

  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #ADJUSTED ARC-SINE METHOD
  for(i in 1:k)
  {
    pAA[i]=y[i]/n1
    qAA[i]=1-pAA[i]
    seAA[i]=cv/sqrt(4*n1)
    LAA[i]=(sin(asin(sqrt(pAA[i]))-seAA[i]))^2
    UAA[i]=(sin(asin(sqrt(pAA[i]))+seAA[i]))^2
    if(LAA[i]<0) LAA[i]=0
    if(UAA[i]>1) UAA[i]=1
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LAA[i] && hp[j] < UAA[i])
      {
        cpAA[i,j]=dbinom(i-1, n,hp[j])
        ctAA[i,j]=1
      }
    }
    cppAA[j]=sum(cpAA[,j])						#Coverage Probability
    if(t1<cppAA[j]&&cppAA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPAA=data.frame(hp,cp=cppAA,method="Adj-ArcSine")
  return(CPAA)
}

##### 4.ADJUSTED LOGIT-WALD - Coverage Probability
gcovpALT<-function(n,alp,h,a,b,t1,t2)
{
  ####INPUT n
  x=0:n
  k=n+1
  y=x+h
  n1=n+(2*h)

  ####INITIALIZATIONS
  pALT=0
  qALT=0
  seALT=0
  lgit=0
  LALT=0
  UALT=0
  s=5000								#Simulation run to generate hypothetical p
  cpALT=matrix(0,k,s)
  ctALT=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppALT=0								#Coverage probabilty
  ctr=0

  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #ADJUSTED LOGIT-WALD METHOD
  for(i in 1:k)
  {
    pALT[i]=y[i]/n1
    qALT[i]=1-pALT[i]
    lgit[i]=log(pALT[i]/qALT[i])
    seALT[i]=sqrt(pALT[i]*qALT[i]*n1)
    LALT[i]=1/(1+exp(-lgit[i]+(cv/seALT[i])))
    UALT[i]=1/(1+exp(-lgit[i]-(cv/seALT[i])))
    if(LALT[i]<0) LALT[i]=0
    if(UALT[i]>1) UALT[i]=1
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LALT[i] && hp[j] < UALT[i])
      {
        cpALT[i,j]=dbinom(i-1, n,hp[j])
        ctALT[i,j]=1
      }
    }
    cppALT[j]=sum(cpALT[,j])				#Coverage Probability
    if(t1<cppALT[j]&&cppALT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

  }
  CPALT=data.frame(hp,cp=cppALT,method="Adj-Logit-Wald")
  return(CPALT)
}

##### 5. ADJUSTED WALD_t - Coverage Probability
gcovpATW<-function(n,alp,h,a,b,t1,t2)
{
  ####INPUT n
  x=0:n
  k=n+1
  y=x+h
  n1=n+(2*h)

  ####INITIALIZATIONS
  pATW=0
  qATW=0
  seATW=0
  LATW=0
  UATW=0
  DOF=0
  cv=0
  s=5000								#Simulation run to generate hypothetical p
  cpATW=matrix(0,k,s)
  ctATW=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppATW=0
  ctr=0

  #MODIFIED_t-WALD METHOD
  for(i in 1:k)
  {
    pATW[i]=y[i]/n1
    qATW[i]=1-pATW[i]
    f1=function(p,n) p*(1-p)/n
    f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
    DOF[i]=2*((f1(pATW[i],n1))^2)/f2(pATW[i],n1)
    cv[i]=qt(1-(alp/2), df=DOF[i])
    seATW[i]=cv[i]*sqrt(f1(pATW[i],n1))
    LATW[i]=pATW[i]-(seATW[i])
    UATW[i]=pATW[i]+(seATW[i])
    if(LATW[i]<0) LATW[i]=0
    if(UATW[i]>1) UATW[i]=1
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LATW[i] && hp[j] < UATW[i])
      {
        cpATW[i,j]=dbinom(i-1, n,hp[j])
        ctATW[i,j]=1
      }
    }
    cppATW[j]=sum(cpATW[,j])						#Coverage Probability
    if(t1<cppATW[j]&&cppATW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPATW=data.frame(hp,cp=cppATW,method="Adj-Wald-T")
  return(CPATW)
}

##### 6.ADJUSTED LIKELIHOOD RATIO - Coverage Probability
gcovpALR<-function(n,alp,h,a,b,t1,t2)
{
  ####INPUT n
  y=0:n
  k=n+1
  y1=y+h
  n1=n+(2*h)
  ####INITIALIZATIONS
  mle=0
  cutoff=0
  LAL=0
  UAL=0
  s=5000								#Simulation run to generate hypothetical p
  cpAL=matrix(0,k,s)
  ctAL=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppAL=0								#Coverage probabilty
  ctr=0

  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #ADJUSTED LIKELIHOOD-RATIO METHOD
  for(i in 1:k)
  {
    likelhd = function(p) dbinom(y1[i],n1,p)
    loglik = function(p) dbinom(y1[i],n1,p,log=TRUE)
    mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
    cutoff[i]=loglik(mle[i])-(cv^2/2)
    loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
    LAL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
    UAL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LAL[i] && hp[j] < UAL[i])
      {
        cpAL[i,j]=dbinom(i-1, n,hp[j])
        ctAL[i,j]=1
      }
    }
    cppAL[j]=sum(cpAL[,j])						#Coverage Probability
    if(t1<cppAL[j]&&cppAL[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPAL=data.frame(hp,cp=cppAL,method="Adj-Likelihood")
  return(CPAL)
}

