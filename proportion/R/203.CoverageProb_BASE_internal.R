#1.WALD
##### 1.WALD-Coverage Probability
gcovpW<-function(n,alp,a,b,t1,t2)
{
  ###INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pW=0
  qW=0
  seW=0
  LW=0
  UW=0
  s=5000								#Simulation run to generate hypothetical p
  cpW=matrix(0,k,s)
  ctW=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppW=0								#Coverage probabilty
  ctr=0
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
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LW[i] && hp[j] < UW[i])
      {
        cpW[i,j]=dbinom(i-1, n,hp[j])
        ctW[i,j]=1
      }
    }
    cppW[j]=sum(cpW[,j])
    if(t1<cppW[j]&&cppW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPW=data.frame(hp,cp=cppW,method="Wald")

  return(CPW)
}

##### 2.SCORE - Coverage Probability
gcovpS<-function(n,alp,a,b,t1,t2)
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
  s=5000								#Simulation run to generate hypothetical p
  cpS=matrix(0,k,s)
  ctS=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppS=0								#Coverage probabilty
  ctr=0

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
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LS[i] && hp[j] < US[i])
      {
        cpS[i,j]=dbinom(i-1, n,hp[j])
        ctS[i,j]=1
      }
    }
    cppS[j]=sum(cpS[,j])						#Coverage Probability
    if(t1<cppS[j]&&cppS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPS=data.frame(hp,cp=cppS,method="Score")
  return(CPS)
}

##### 3.ARC SINE - Coverage Probability
gcovpA<-function(n,alp,a,b,t1,t2)
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
  s=5000								#Simulation run to generate hypothetical p
  cpA=matrix(0,k,s)
  ctA=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppA=0								#Coverage probabilty
  ctr=0

  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
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
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LA[i] && hp[j] < UA[i])
      {
        cpA[i,j]=dbinom(i-1, n,hp[j])
        ctA[i,j]=1
      }
    }
    cppA[j]=sum(cpA[,j])						#Coverage Probability
    if(t1<cppA[j]&&cppA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPA=data.frame(hp,cp=cppA,method="ArcSine")
  return(CPA)
}

##### 4.LOGIT-WALD - Coverage Probability
gcovpLT<-function(n,alp,a,b,t1,t2)
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
  s=5000								#Simulation run to generate hypothetical p
  cpLT=matrix(0,k,s)
  ctLT=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppLT=0								#Coverage probabilty
  ctr=0

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
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LLT[i] && hp[j] < ULT[i])
      {
        cpLT[i,j]=dbinom(i-1, n,hp[j])
        ctLT[i,j]=1
      }
    }
    cppLT[j]=sum(cpLT[,j])				#Coverage Probability
    if(t1<cppLT[j]&&cppLT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

  }
  CPLT=data.frame(hp,cp=cppLT,method="Logit-Wald")
  return(CPLT)
}

##### 5. WALD_t - Coverage Probability
gcovpTW<-function(n,alp,a,b,t1,t2)
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
  s=5000								#Simulation run to generate hypothetical p
  cpTW=matrix(0,k,s)
  ctTW=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppTW=0
  ctr=0

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
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LTW[i] && hp[j] < UTW[i])
      {
        cpTW[i,j]=dbinom(i-1, n,hp[j])
        ctTW[i,j]=1
      }
    }
    cppTW[j]=sum(cpTW[,j])						#Coverage Probability
    if(t1<cppTW[j]&&cppTW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPTW=data.frame(hp,cp=cppTW,method="Wald-T")
  return(CPTW)
}

##### 6.LIKELIHOOD RATIO - Coverage Probability
gcovpL<-function(n,alp,a,b,t1,t2)
{
  ####INPUT n
  y=0:n
  k=n+1
  ####INITIALIZATIONS
  mle=0
  cutoff=0
  LL=0
  UL=0
  s=5000								#Simulation run to generate hypothetical p
  cpL=matrix(0,k,s)
  ctL=matrix(0,k,s)							#Cover Pbty quantity in sum
  cppL=0								#Coverage probabilty
  ctr=0

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
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LL[i] && hp[j] < UL[i])
      {
        cpL[i,j]=dbinom(i-1, n,hp[j])
        ctL[i,j]=1
      }
    }
    cppL[j]=sum(cpL[,j])						#Coverage Probability
    if(t1<cppL[j]&&cppL[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPL=data.frame(hp,cp=cppL,method="Likelihood")
  return(CPL)
}

##### 7.Exact - Coverage Probability
gcovpEX=function(n,alp,e,a,b,t1,t2)
{
  nvar=length(e)

  res <- data.frame()

  for(i in 1:nvar)
  {
    lu=gintcovpEX202(n,alp,e[i],a,b,t1,t2)
    res <- rbind(res,lu)
  }
  return(res)
}
gintcovpEX202=function(n,alp,e,a,b,t1,t2)
{

  x=0:n
  k=n+1
  LEX=0
  UEX=0
  s=5000					#Simulation run to generate hypothetical p
  cpEX=matrix(0,k,s)
  ctEX=matrix(0,k,s)			#Cover Pbty quantity in sum
  cppEX=0
  ctr=0
  #EXACT METHOD
  LEX[1]=0
  UEX[1]= 1-((alp/(2*e))^(1/n))
  LEX[k]=(alp/(2*e))^(1/n)
  UEX[k]=1

  for(i in 1:(k-2))
  {
    LEX[i+1]=exlim202l(x[i+1],n,alp,e)
    UEX[i+1]=exlim202u(x[i+1],n,alp,e)
  }
  ####COVERAGE PROBABILITIES
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      if(hp[j] > LEX[i] && hp[j] < UEX[i])
      {
        cpEX[i,j]=dbinom(i-1, n,hp[j])
        ctEX[i,j]=1
      }
    }
    cppEX[j]=sum(cpEX[,j])						#Coverage Probability
    if(t1<cppEX[j]&&cppEX[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
  }
  CPEX=data.frame(hp,cpp=cppEX)
  mcpEX=mean(cppEX)							#Mean Cov Prob
  micpEX=min(cppEX)							#Min Cov Prob
  return(data.frame(CPEX,mcpEX,micpEX,e))
}
#####TO FIND LOWER LIMITS
exlim202l=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
  LEX= uniroot(f1,c(0,1))$root
  return(LEX)
}
#####TO FIND UPPER LIMITS
exlim202u=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f2  = function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
  UEX = uniroot(f2,c(0,1))$root
  return(UEX)
}
