##### 1.CC-WALD Expected Length for a given n and alpha level
gexplCWD<-function(n,alp,c,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pCW=0
  qCW=0
  seCW=0
  LCW=0
  UCW=0
  s=5000
  LECW=0 								#LENGTH OF INTERVAL

  ewiCW=matrix(0,k,s)						#Expected length quantity in sum
  ewCW=0									#Expected Length
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #WALD METHOD
  for(i in 1:k)
  {
    pCW[i]=x[i]/n
    qCW[i]=1-pCW[i]
    seCW[i]=sqrt(pCW[i]*qCW[i]/n)
    LCW[i]=pCW[i]-((cv*seCW[i])+c)
    UCW[i]=pCW[i]+((cv*seCW[i])+c)
    if(LCW[i]<0) LCW[i]=0
    if(UCW[i]>1) UCW[i]=1
    LECW[i]=UCW[i]-LCW[i]
  }
  #sumLECW=sum(LECW)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiCW[i,j]=LECW[i]*dbinom(i-1, n,hp[j])
    }
    ewCW[j]=sum(ewiCW[,j])						#Expected Length
  }
  ELCW=data.frame(hp,ew=ewCW,method="CC-Wald")
  return(ELCW)
}

##### 2.CC SCORE - Expected Length for a given n and alpha level
gexplCSC<-function(n,alp,c,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pCS=0
  qCS=0
  seCS_L=0
  seCS_U=0
  LCS=0
  UCS=0
  s=5000
  LECS=0 								#LENGTH OF INTERVAL

  ewiCS=matrix(0,k,s)						#Expected length quantity in sum
  ewCS=0									#Expected Length
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  cv1=(cv^2)/(2*n)
  cv2= cv/(2*n)

  #SCORE (WILSON) METHOD
  for(i in 1:k)
  {
    pCS[i]=x[i]/n
    qCS[i]=1-pCS[i]
    seCS_L[i]=sqrt((cv^2)-(4*n*(c+c^2))+(4*n*pCS[i]*(1-pCS[i]+(2*c))))	#Sq. root term of LL
    seCS_U[i]=sqrt((cv^2)+(4*n*(c-c^2))+(4*n*pCS[i]*(1-pCS[i]-(2*c))))	#Sq. root term of LL
    LCS[i]=(n/(n+(cv)^2))*((pCS[i]-c+cv1)-(cv2*seCS_L[i]))
    UCS[i]=(n/(n+(cv)^2))*((pCS[i]+c+cv1)+(cv2*seCS_U[i]))
    if(LCS[i]<0) LCS[i]=0
    if(UCS[i]>1) UCS[i]=1
    LECS[i]=UCS[i]-LCS[i]
  }
  #sumLECS=sum(LECS)

  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiCS[i,j]=LECS[i]*dbinom(i-1, n,hp[j])
    }
    ewCS[j]=sum(ewiCS[,j])						#Expected Length
  }
  ELCS=data.frame(hp,ew=ewCS,method="CC-Wilson")
  return(ELCS)
}

##### 3.CC ARC SINE - Expected Length for a given n and alpha level
gexplCAS<-function(n,alp,c,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pCA=0
  qCA=0
  seCA=0
  LCA=0
  UCA=0
  s=5000
  LECA=0 								#LENGTH OF INTERVAL

  ewiCA=matrix(0,k,s)						#Expected length quantity in sum
  ewCA=0									#Expected Length
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #ARC-SINE METHOD
  for(i in 1:k)
  {
    pCA[i]=x[i]/n
    qCA[i]=1-pCA[i]
    seCA[i]=cv/sqrt(4*n)
    LCA[i]=(sin(asin(sqrt(pCA[i]))-seCA[i]-c))^2
    UCA[i]=(sin(asin(sqrt(pCA[i]))+seCA[i]+c))^2
    if(LCA[i]<0) LCA[i]=0
    if(UCA[i]>1) UCA[i]=1
    LECA[i]=UCA[i]-LCA[i]
  }
  #sumLECA=sum(LECA)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiCA[i,j]=LECA[i]*dbinom(i-1, n,hp[j])
    }
    ewCA[j]=sum(ewiCA[,j])						#Expected Length
  }
  ELCA=data.frame(hp,ew=ewCA,method="CC-ArcSine")
  return(ELCA)
}

##### 4.CC LOGIT-WALD - Expected Length for a given n and alpha level
gexplCLT<-function(n,alp,c,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  #Expected Length
  #INITIALIZATIONS
  pCLT=0
  qCLT=0
  seCLT=0
  lgit=0
  LCLT=0
  UCLT=0
  s=5000
  LECLT=0 								#LENGTH OF INTERVAL

  ewiCLT=matrix(0,k,s)						#Expected length quantity in sum
  ewCLT=0									#Expected Length
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #LOGIT-WALD METHOD
  pCLT[1]=0
  qCLT[1]=1
  LCLT[1] = 0
  UCLT[1] = 1-((alp/2)^(1/n))

  pCLT[k]=1
  qCLT[k]=0
  LCLT[k]= (alp/2)^(1/n)
  UCLT[k]=1

  lgiti=function(t) exp(t)/(1+exp(t))	#LOGIT INVERSE
  for(j in 1:(k-2))
  {
    pCLT[j+1]=x[j+1]/n
    qCLT[j+1]=1-pCLT[j+1]
    lgit[j+1]=log(pCLT[j+1]/qCLT[j+1])
    seCLT[j+1]=sqrt(pCLT[j+1]*qCLT[j+1]*n)
    LCLT[j+1]=lgiti(lgit[j+1]-(cv/seCLT[j+1])-c)
    UCLT[j+1]=lgiti(lgit[j+1]+(cv/seCLT[j+1])+c)
  }
  for(i in 1:k)
  {
    if(LCLT[i]<0) LCLT[i]=0
    if(UCLT[i]>1) UCLT[i]=1
    LECLT[i]=UCLT[i]-LCLT[i]
  }
  #sumLECLT=sum(LECLT)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiCLT[i,j]=LECLT[i]*dbinom(i-1, n,hp[j])
    }
    ewCLT[j]=sum(ewiCLT[,j])						#Expected Length
  }
  ELCLT=data.frame(hp,ew=ewCLT,method="CC-Logit-Wald")
  return(ELCLT)
}

##### 5.CC t-WALD_CC - Expected Length for a given n and alpha level
gexplCTW<-function(n,alp,c,a,b)
{
  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pCTW=0
  qCTW=0
  seCTW=0
  LCTW=0
  UCTW=0
  DOF=0
  cv=0
  s=5000
  LECTW=0 								#LENGTH OF INTERVAL

  ewiCTW=matrix(0,k,s)						#Expected length quantity in sum
  ewCTW=0									#Expected Length
  #MODIFIED_t-WALD METHOD
  for(i in 1:k)
  {
    if(x[i]==0||x[i]==n)
    {
      pCTW[i]=(x[i]+2)/(n+4)
      qCTW[i]=1-pCTW[i]
    }else
    {
      pCTW[i]=x[i]/n
      qCTW[i]=1-pCTW[i]
    }
    f1=function(p,n) p*(1-p)/n
    f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
    DOF[i]=2*((f1(pCTW[i],n))^2)/f2(pCTW[i],n)
    cv[i]=qt(1-(alp/2), df=DOF[i])
    seCTW[i]=cv[i]*sqrt(f1(pCTW[i],n))
    LCTW[i]=pCTW[i]-(seCTW[i]+c)
    UCTW[i]=pCTW[i]+(seCTW[i]+c)
    if(LCTW[i]<0) LCTW[i]=0
    if(UCTW[i]>1) UCTW[i]=1
    LECTW[i]=UCTW[i]-LCTW[i]
  }
  #sumLECTW=sum(LECTW)
  ####Expected Length
  hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
  for (j in 1:s)
  {
    for(i in 1:k)
    {
      ewiCTW[i,j]=LECTW[i]*dbinom(i-1, n,hp[j])
    }
    ewCTW[j]=sum(ewiCTW[,j])						#Expected Length
  }
  ELCTW=data.frame(hp,ew=ewCTW,method="CC-Wald-T")
  return(ELCTW)
}

