 #rm(list=ls(all=TRUE))

BiGQD.mle=function(X,time,mesh=10,theta,control=NULL,method='Nelder-Mead',RK.order=4,exclude=NULL,Tag=NA,Dtype='Saddlepoint',rtf= runif(2,-1,1),wrt=FALSE,print.output=TRUE)
{
  rtf=0*rtf+1
  solver   =function(Xs, Xt, theta, N , delt , N2, tt  , P , alpha, lower , upper, tro  ){}
  rm(list =c('solver')) 
  check_for_model=function()
  {
    txt=''
    namess=c('a00','a10','a20','a01','a02','a11',
           'b00','b10','b20','b01','b02','b11',
           'c00','c10','c20','c01','c02','c11',
           'd00','d10','d20','d01','d02','d11',
           'e00','e10','e20','e01','e02','e11',
           'f00','f10','f20','f01','f02','f11')
    func.list=rep(0,length(namess))
    obs=objects(pos=1)
    for(i in 1:length(namess))
    {
    if(sum(obs==namess[i])){func.list[i]=1}
    }

    check=F
    if(sum(func.list)==0)
   {
   txt='
  --------------------------------------------------------------------------------
  No model has been defined yet! Try for example:
  --------------------------------------------------------------------------------
  GQD.remove()
  a00=function(t){theta[1]*theta[2]}
  a10=function(t){-theta[1]}
  c00=function(t){theta[3]*theta[3]}

  b00=function(t){theta[4]*theta[5]}
  b01=function(t){-theta[4]}
  f00=function(t){theta[6]*theta[6]}

  model=BiGQD.mle(X,time,10,theta =rep(0.1,6))
  --------------------------------------------------------------------------------
   '
   check=T
   }
   if((sum(func.list)>0)&&(sum(func.list[-c(1:12)])==0))
   {
   txt='
  --------------------------------------------------------------------------------
  At least two diffusion coefficients have to be defined! Try for example:
  --------------------------------------------------------------------------------
  GQD.remove()
  c00=function(t){theta[1]*theta[1]}
  f00=function(t){theta[2]*theta[2]}
  model=BiGQD.mle(X,time,10,theta =rep(0.1,2))
  --------------------------------------------------------------------------------
   '
     check=T
    }
    return(list(check=check,txt=txt))
  }
  check_for=check_for_model()
  if(check_for[[1]]){stop(check_for[[2]])}
  theta = theta+runif(length(theta),0.001,0.002)*sign(theta)
  T.seq=time
  Dtypes =c('Saddlepoint','Edgeworth','Normal')
  Dindex = which(Dtypes==Dtype)

  # A power function used to solve the syntax gap issue between R and C++
  pow=function(x,p)
  {
    x^p
  }
  prod=function(a,b){a*b}

  # Separate the time series into X and Y components
  X1=X[,1]
  X2=X[,2]
  nnn=length(X1)

   b1 = '\n==============================================================================\n'
   b2 = '==============================================================================\n'
  warn=c(
    '1.  Missing input: Argument {X} is missing.\n'
    ,'2.  Missing input: Argument {time} is missing.\n'
    ,'3.  Missing input: Argument {theta} is missing.\n'
    ,'4.  Missing input: Argument {sds} is missing.\n'
    ,'5.  Input type: Argument {X} must be of type matrix!.\n'
    ,'6.  Input type: Argument {time} must be of type vector!.\n'
    ,'7.  Input: Less starting parameters than model parameters.\n'
    ,'8.  Input: More starting parameters than model parameters.\n'
    ,'9.  Input: length(X) must be > 10.\n'
    ,'10. Input: length(time) must be > 10.\n'
    ,'11. Input: length(lower)!=1.\n'
    ,'12. Input: length(upper)!=1.\n'
    ,'13. Input: length(P)!=1.\n'
    ,'14. Input: length(mesh)!=1.\n'
    ,'15. Input: length(alpha)!=1.\n'
    ,'16. Input: length(Trunc)!=1.\n'
    ,'17. Input: length(RK.order)!=1.\n'
    ,'18. Density: Dtype has to be one of Saddlepoint, Edgeworth or Normal.\n'
    ,'19. Density: Range [lower,upper] must be strictly positive for Dtype Gamma or InvGamma.\n'
    ,'20. Density: Dtype cannot be Beta for observations not in (0,1).\n'
    ,'21. Density: Argument {upper} must be > {lower}.\n'
    ,'22. Density: P must be >= 10.\n'
    ,'23. Density: Trunc[2] must be <= Trunc[1].\n'
    ,'24. ODEs : Large max(diff(time))/mesh may result in poor approximations. Try larger mesh.\n'
    ,'25. ODEs : max(diff(time))/mesh must be <1.\n'
    ,'26. ODEs : Runge Kutta scheme must be of order 4 or 10.\n'
    ,'27. ODEs : Argument {mesh} must be >= 5.\n'
    ,'28. Input: length(X)!=length(time).\n'
    ,'29. MCMC : Argument {burns} must be < {N.updates}.\n'
    ,'30. MCMC : Argument {N.updates} must be > 2.\n'
    ,'31. MCMC : length(theta)!=length(sds).\n'
    ,'32. Model: There has to be at least one model coefficient.\n'
    ,'33. Input: length(N.updates)!=1.\n'
    ,'34. Input: length(burns)!=1.\n'
    ,'35. Prior: priors(theta) must return a single value.\n'
    ,'36. Input: NAs not allowed.\n'
    ,'37. Input: length(Dtype)!=1.\n'
    ,'38. Input: NAs not allowed.\n'
    ,'39. Input: Time series contains values of small magnitude.\n    This may result in numerical instabilities.\n    It may be advisable to scale the data by a constant factor.\n'
  )

   warntrue = rep(F,40)
    check.thetas = function(theta,tt)
  {
    t=tt
    theta = theta+runif(length(theta),0.001,0.002)*sign(theta)
    namess=c('a00','a10','a20','a01','a02','a11',
             'b00','b10','b20','b01','b02','b11',
             'c00','c10','c20','c01','c02','c11',
             'd00','d10','d20','d01','d02','d11',
             'e00','e10','e20','e01','e02','e11',
             'f00','f10','f20','f01','f02','f11')
    func.list=rep(0,length(namess))
    obs=objects(pos=1)
    for(i in 1:length(namess))
    {
      if(sum(obs==namess[i])){func.list[i]=1}
    }
    pers.represented = rep(0,length(theta))
    for(i in which(func.list==1))
    {
      for(j in 1:length(theta))
      {
        dresult1=eval(body(namess[i]))
        theta[j] = theta[j]+runif(1,0.1,0.2)
        dresult2=eval(body(namess[i]))
        dff = abs(dresult1-dresult2)
        if(any(round(dff,6)!=0)){pers.represented[j]=pers.represented[j]+1}
      }
    }
    return(pers.represented)
  }
  check.thetas2 = function(theta)
  {
    namess=c('a00','a10','a20','a01','a02','a11',
             'b00','b10','b20','b01','b02','b11',
             'c00','c10','c20','c01','c02','c11',
             'd00','d10','d20','d01','d02','d11',
             'e00','e10','e20','e01','e02','e11',
             'f00','f10','f20','f01','f02','f11')
    func.list=rep(0,length(namess))
    obs=objects(pos=1)
    for(i in 1:length(namess))
    {
      if(sum(obs==namess[i])){func.list[i]=1}
    }

    l=0
    for(k in which(func.list==1))
    {
      str=body(namess[k])[2]
      for(i in 1:length(theta))
      {
        for(j in 1:20)
        {
          str=sub(paste0('theta\\[',i,'\\]'),'clear',str)
        }
      }
      l=l+length(grep('theta',str))
      l
    }
    return(l)
  }

   # Check missing values first:
   if(missing(X))                                                {warntrue[1]=T}
   if(missing(time))                                             {warntrue[2]=T}
   if(missing(theta))                                            {warntrue[3]=T}
   #if(missing(sds))                                              {warntrue[4]=T}
   if(!is.matrix(X))                                             {warntrue[5]=T}
   if(!is.vector(time))                                          {warntrue[6]=T}
   # Check model parameters:
   if(check.thetas2(theta)!=0)                                   {warntrue[7]=T}
   if(!warntrue[7]){if(any(check.thetas(theta,T.seq)==0))        {warntrue[8]=T}}

   # Check input length:
   if(dim(X)[1]<10)                                              {warntrue[9]=T}
   if(length(time)<10)                                          {warntrue[10]=T}
   #if(length(lower)>1)                                          {warntrue[11]=T}
   #if(length(upper)>1)                                          {warntrue[12]=T}
   #if(length(P)!=1)                                             {warntrue[13]=T}
   if(length(mesh)!=1)                                        {warntrue[14]=T}
   #if(length(alpha)!=1)                                         {warntrue[15]=T}
   #if(length(Trunc)!=2)                                         {warntrue[16]=T}
   if(length(RK.order)!=1)                                      {warntrue[17]=T}
   #if(length(N.updates)!=1)                                     {warntrue[33]=T}
   #if(length(burns)!=1)                                         {warntrue[34]=T}
   if(length(Dtype)!=1)                                         {warntrue[37]=T}


   # Check density approx parameters:
   if(sum(Dindex)==0)                                          {warntrue[18] =T}
   #if(!warntrue[18])
   #{
   # if((Dindex==3)|(Dindex==4)){if(lower[1]<=0)                {warntrue[19] =T}}
   # if(Dindex==5){if(any(X<=0)|any(X>=1))                      {warntrue[20] =T}}
   #}
   #if(!any(warntrue[c(11,12)])){if(upper<=lower)                {warntrue[21] =T}}
   #if(!warntrue[13]){if(P<10)                                   {warntrue[22] =T}}
   #if(!warntrue[16]){if(Trunc[2]>Trunc[1])                      {warntrue[23] =T}}

   #  Miscelaneous checks:
   excl=0
   if(is.null(exclude)){excl=length(T.seq)-1+200}
   if(!is.null(exclude)){excl=exclude}
   test.this =max(diff(T.seq)[-excl])/mesh
   if(test.this>0.1)                                            {warntrue[24]=T}
   if(test.this>=1)                                             {warntrue[25]=T}
   if(!warntrue[17]){if(!((RK.order==4)|(RK.order==10)))        {warntrue[26]=T}}
   if(!warntrue[14]){if(mesh<5)                               {warntrue[27]=T}}
   if(dim(X)[1]!=length(time))                                  {warntrue[28]=T}
   #if(!any(warntrue[c(33,34)])){if(burns>N.updates)             {warntrue[29]=T}}
   #if(!warntrue[33]){if(N.updates<2)                            {warntrue[30]=T}}
   #if(length(theta)!=length(sds))                               {warntrue[31]=T}
   if(any(is.na(X))||any(is.na(time)))                          {warntrue[36]=T}

   # Print output:
   if(any(warntrue))
   {
      prnt = b1
      for(i in which(warntrue))
      {
         prnt = paste0(prnt,warn[i])
      }
      prnt = paste0(prnt,b2)
      stop(prnt)
   }
     if(any(X<10^-2)){warntrue[39]=T}
   # Print warnings:
   if(any(warntrue))
   {
      prnt = b1
      for(i in which(warntrue))
      {
         prnt = paste0(prnt,warn[i])
      }
      prnt = paste0(prnt,b2)
      warning(prnt)
   }
    ############################################################################
    ############################################################################
    ############################################################################


   #==============================================================================
   #                           Data resolution Module
   #==============================================================================
   #t=seq(0,100,by=1/100)
   homo=T
   homo.res=T
   delt=(diff(T.seq)/mesh)[1]
   t=T.seq

    if(is.null(exclude))
    {
       if(sum(round(diff(T.seq)-diff(T.seq)[1],10)==0)!=length(T.seq)-1){homo.res=F}
    }
    if(!is.null(exclude))
    {
       if(sum(round(diff(T.seq)[-excl]-c(diff(T.seq)[-excl])[1],10)==0)!=length(T.seq)-1-length(excl)){homo.res=F}
    }

   #==============================================================================
   #                           Solution TYPE Module
   #==============================================================================
    namess=c('a00','a10','a20','a01','a02','a11',
             'b00','b10','b20','b01','b02','b11',
             'c00','c10','c20','c01','c02','c11',
             'd00','d10','d20','d01','d02','d11',
             'e00','e10','e20','e01','e02','e11',
             'f00','f10','f20','f01','f02','f11')
    func.list=rep(0,36)
    obs=objects(pos=1)
    for(i in 1:36)
    {
      if(sum(obs==namess[i])){func.list[i]=1}
    }

    set1.1=c(3,14,15,20,21)                          # X Nonlinear
    set1.2=c(11,28,29,34,35)                         # Y Nonlinear
    set4=c(3,5,6,9,11,12,14:18,20:24,26:30,32:36)    # Nonlinear terms
    set5=c(4:6,8,9,12,16:18,22:24,26,27,30,32,33,36) # Dependant terms

    state3.0=(     any(func.list[set1.1]==1)     &     any(func.list[set1.2]==1)     &(sum(func.list[set5])==0))  # Any nonlinear in BOTH     and NO interaction term
    state3.1=((sum(any(func.list[set1.1]==1))==0)&(    any(func.list[set1.2]==1)    )&(sum(func.list[set5])==0))  # Any nonlinear Y ,Normal X and NO interaction term
    state3.2=(     any(func.list[set1.1]==1)     &(sum(any(func.list[set1.2]==1))==0)&(sum(func.list[set5])==0))  # Any nonlinear X ,Normal Y and NO interaction term
    state4  =(any(func.list[set4]==1)&any(func.list[set5]==1))                                                    # Any nonlinear term AND interaction terms
    state1  =!(state3.0|state3.1|state3.2|state4)
    sol.state=c(state1,state3.0,state3.1,state3.2,state4)
    sol.state=c('2nd Ord. Truncation, Bivariate Normal','4th Ord. Truncation, Independent Saddlepoints ','2nd & 4th Ord. Truncation (Indep)','4th & 2nd Ord. Trunc. (Indep)',"4th Ord. Truncation, Bivariate-Saddlepoint")[which(sol.state==T)]

   #==============================================================================
   #                    Check the syntax of each coefficient
   #==============================================================================
    strcheck=function(str)
    {
      k=0
      if(length(grep('/',str))!=0)
      {
       warning('C++ uses integer division when denominators are of type int. Use decimal values if possible (e.g. 0.5 i.s.o. 1/2).',call. = F)
       k=k+1
      }
      if(length(grep('\\^',str))!=0)
      {
       stop('^-Operator not overloaded in C++: Use pow(x,p) instead (e.g. pow(x,2) i.s.o. x^2).',call. = F)
       k=k+1
      }
      return(k)
    }
    for(i in which(func.list==1))
    {
      strcheck(body(namess[i])[2])
    }

    ############################################################################
    ############################################################################
    ############################################################################

    tro=c(5,8,8,8,14)[which(c(state1,state3.0,state3.1,state3.2,state4)==T)]

if(tro==5)
{
    txtA='
    #include <RcppArmadillo.h>
    #include <math.h>
    #include <Rcpp.h>
    #define pi           3.14159265358979323846  /* pi */
    using namespace arma;
    using namespace Rcpp;
    using namespace R;

    // [[Rcpp::depends("RcppArmadillo")]]
    // [[Rcpp::export]]
    vec prod(vec a,vec b)
    {
        return(a%b);
    }
    mat f(mat a,vec theta,vec t,int N2)
    {
        mat atemp(N2,5);'

secmom = 3
}
if(tro==8)
{
    txtA='
    #include <RcppArmadillo.h>
    #include <math.h>
    #include <Rcpp.h>
    #define pi           3.14159265358979323846  /* pi */
    using namespace arma;
    using namespace Rcpp;
    using namespace R;

    // [[Rcpp::depends("RcppArmadillo")]]
    // [[Rcpp::export]]
    vec prod(vec a,vec b)
    {
        return(a%b);
    }
    mat f(mat a,vec theta,vec t,int N2)
    {
        mat atemp(N2,8);'

secmom = 5
}

if(tro==14)
{
    txtA='
    #include <RcppArmadillo.h>
    #include <math.h>
    #include <Rcpp.h>
    #define pi           3.14159265358979323846  /* pi */
    using namespace arma;
    using namespace Rcpp;
    using namespace R;

    // [[Rcpp::depends("RcppArmadillo")]]
    // [[Rcpp::export]]
    vec prod(vec a,vec b)
    {
        return(a%b);
    }
    mat f(mat a,vec theta,vec t,int N2)
    {
        mat atemp(N2,14);'
secmom = 5
}

if(RK.order==4)
{
  if(homo.res)
  {
txtB= '
    return atemp;
}

// [[Rcpp::export]]
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N,double delt,int N2,vec tt,mat starts,int tro,int secmom)
{
    mat resss(N2,3);
    mat x0(N2,tro);
    mat xa(N2,tro);
    mat xe(N2,tro);
    mat fx1(N2,tro);
    mat fx2(N2,tro);
    mat fx3(N2,tro);
    mat fx4(N2,tro);
    mat fx5(N2,tro);
    mat fx6(N2,tro);
    double whch =0;
    x0.fill(0);
    x0.col(0)=Xs;
    x0.col(secmom-1)=Ys;
    vec d=tt;
    for (int i = 1; i < N; i++)
   {
    fx1 = f(x0,theta,d,N2)*delt;
    fx2 = f(x0+0.25*fx1,theta,d+0.25*delt,N2)*delt;
    fx3 = f(x0+0.09375*fx1+0.28125*fx2,theta,d+0.375*delt,N2)*delt;
    fx4 = f(x0+0.879381*fx1-3.277196*fx2+ 3.320892*fx3,theta,d+0.9230769*delt,N2)*delt;
    fx5 = f(x0+2.032407*fx1-8*fx2+7.173489*fx3-0.2058967*fx4,theta,d+delt,N2)*delt;
    fx6 = f(x0-0.2962963*fx1+2*fx2-1.381676*fx3+0.4529727*fx4-0.275*fx5,theta,d+0.5*delt,N2)*delt;
    xa  =  x0+0.1185185*fx1+0.5189864*fx3+0.5061315*fx4-0.18*fx5+0.03636364*fx6;
    x0  =  x0+0.1157407*fx1+0.5489279*fx3+0.5353314*fx4-0.2*fx5;
    xe  =  abs(x0.col(1)-xa.col(1));
    if(xe.max()>whch)
    {
      whch = xe.max();
    }
    d=d+delt;
   }
    '
  }
  if(!homo.res)
  {
txtB= '
    return atemp;
}

// [[Rcpp::export]]
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N, mat delt,int N2,vec d,mat starts,int tro,int secmom)
{
    mat resss(N2,3);
    mat x0(N2,tro);
    mat xa(N2,tro);
    mat xe(N2,tro);
    mat fx1(N2,tro);
    mat fx2(N2,tro);
    mat fx3(N2,tro);
    mat fx4(N2,tro);
    mat fx5(N2,tro);
    mat fx6(N2,tro);
    double whch =0;
    x0.fill(0);
    x0.col(0)=Xs;
    x0.col(secmom-1)=Ys;
    for (int i = 1; i < N; i++)
   {
    fx1 = f(x0,theta,d,N2)%delt;
    fx2 = f(x0+0.25*fx1,theta,d+0.25*delt.col(0),N2)%delt;
    fx3 = f(x0+0.09375*fx1+0.28125*fx2,theta,d+0.375*delt.col(0),N2)%delt;
    fx4 = f(x0+0.879381*fx1-3.277196*fx2+ 3.320892*fx3,theta,d+0.9230769*delt.col(0),N2)%delt;
    fx5 = f(x0+2.032407*fx1-8*fx2+7.173489*fx3-0.2058967*fx4,theta,d+delt.col(0),N2)%delt;
    fx6 = f(x0-0.2962963*fx1+2*fx2-1.381676*fx3+0.4529727*fx4-0.275*fx5,theta,d+0.5*delt.col(0),N2)%delt;
    xa  =  x0+0.1185185*fx1+0.5189864*fx3+0.5061315*fx4-0.18*fx5+0.03636364*fx6;
    x0  =  x0+0.1157407*fx1+0.5489279*fx3+0.5353314*fx4-0.2*fx5;
    xe  =  abs(x0.col(1)-xa.col(1));
    if(xe.max()>whch)
    {
      whch = xe.max();
    }
    d=d+delt.col(0);
   }
    '
   }
}
if(RK.order==10)
{
  if(homo.res)
  {
txtB= '
  return atemp;
}

// [[Rcpp::export]]
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N,double delt,int N2,vec tt,mat starts,int tro,int secmom)
{
  mat resss(N2,3);
  mat fx0(N2,tro);
  mat fx1(N2,tro);
  mat fx2(N2,tro);
  mat fx3(N2,tro);
  mat fx4(N2,tro);
  mat fx5(N2,tro);
  mat fx6(N2,tro);
  mat fx7(N2,tro);
  mat fx8(N2,tro);
  mat fx9(N2,tro);
  mat fx10(N2,tro);
  mat fx11(N2,tro);
  mat fx12(N2,tro);
  mat fx13(N2,tro);
  mat fx14(N2,tro);
  mat fx15(N2,tro);
  mat fx16(N2,tro);
  mat x0(N2,tro);
  mat x1(N2,tro);
  mat x2(N2,tro);
  mat x3(N2,tro);
  mat x4(N2,tro);
  mat x5(N2,tro);
  mat x6(N2,tro);
  mat x7(N2,tro);
  mat x8(N2,tro);
  mat x9(N2,tro);
  mat x10(N2,tro);
  mat x11(N2,tro);
  mat x12(N2,tro);
  mat x13(N2,tro);
  mat x14(N2,tro);
  mat x15(N2,tro);
  mat x16(N2,tro);
  x0.fill(0);
  x0.col(0)=Xs;
  x0.col(secmom-1)=Ys;
  vec d=tt;
  for (int i = 1; i < N; i++)
{

  fx0=f(x0,theta,d,N2);
  x1=x0+delt*(0.1*fx0);
  fx1=f(x1,theta,d+0.100000000000000000000000000000000000000000000000000000000000*delt,N2);
  x2=x0+delt*(-0.915176561375291*fx0+1.45453440217827*fx1);
  fx2=f(x2,theta,d+0.539357840802981787532485197881302436857273449701009015505500*delt,N2);
  x3=x0+delt*(0.202259190301118*fx0+0.606777570903354*fx2);
  fx3=f(x3,theta,d+0.809036761204472681298727796821953655285910174551513523258250*delt,N2);
  x4=x0+delt*(0.184024714708644*fx0+0.197966831227192*fx2-0.0729547847313633*fx3);
  fx4=f(x4,theta,d+0.309036761204472681298727796821953655285910174551513523258250*delt,N2);
  x5=x0+delt*(0.0879007340206681*fx0+0.410459702520261*fx3+0.482713753678866*fx4);
  fx5=f(x5,theta,d+0.981074190219795268254879548310562080489056746118724882027805*delt,N2);
  x6=x0+delt*(0.085970050490246*fx0+0.330885963040722*fx3+0.48966295730945*fx4-0.0731856375070851*fx5);
  fx6=f(x6,theta,d+0.833333333333333333333333333333333333333333333333333333333333*delt,N2);
  x7=x0+delt*(0.120930449125334*fx0+0.260124675758296*fx4+0.0325402621549091*fx5-0.0595780211817361*fx6);
  fx7=f(x7,theta,d+0.354017365856802376329264185948796742115824053807373968324184*delt,N2);
  x8=x0+delt*(0.110854379580391*fx0-0.0605761488255006*fx5+0.321763705601778*fx6+0.510485725608063*fx7);
  fx8=f(x8,theta,d+0.882527661964732346425501486979669075182867844268052119663791*delt,N2);
  x9=x0+delt*(0.112054414752879*fx0-0.144942775902866*fx5-0.333269719096257*fx6+0.49926922955688*fx7+0.509504608929686*fx8);
  fx9=f(x9,theta,d+0.642615758240322548157075497020439535959501736363212695909875*delt,N2);
  x10=x0+delt*(0.113976783964186*fx0-0.0768813364203357*fx5+0.239527360324391*fx6+0.397774662368095*fx7+0.0107558956873607*fx8-0.327769124164019*fx9);
  fx10=f(x10,theta,d+0.357384241759677451842924502979560464040498263636787304090125*delt,N2);
  x11=x0+delt*(0.0798314528280196*fx0-0.0520329686800603*fx5-0.0576954146168549*fx6+0.194781915712104*fx7+0.145384923188325*fx8-0.0782942710351671*fx9-0.114503299361099*fx10);
  fx11=f(x11,theta,d+0.117472338035267653574498513020330924817132155731947880336209*delt,N2);
  x12=x0+delt*(0.985115610164857*fx0+0.330885963040722*fx3+0.48966295730945*fx4-1.37896486574844*fx5-0.861164195027636*fx6+5.78428813637537*fx7+3.28807761985104*fx8-2.38633905093136*fx9-3.25479342483644*fx10-2.16343541686423*fx11);
  fx12=f(x12,theta,d+0.833333333333333333333333333333333333333333333333333333333333*delt,N2);
  x13=x0+delt*(0.895080295771633*fx0+0.197966831227192*fx2-0.0729547847313633*fx3-0.851236239662008*fx5+0.398320112318533*fx6+3.63937263181036*fx7+1.5482287703983*fx8-2.12221714704054*fx9-1.58350398545326*fx10-1.71561608285936*fx11-0.0244036405750127*fx12);
  fx13=f(x13,theta,d+0.309036761204472681298727796821953655285910174551513523258250*delt,N2);
  x14=x0+delt*(-0.915176561375291*fx0+1.45453440217827*fx1+0*fx2+0*fx3-0.777333643644968*fx4+0*fx5-0.0910895662155176*fx6+0.0910895662155176*fx12+0.777333643644968*fx13);
  fx14=f(x14,theta,d+0.539357840802981787532485197881302436857273449701009015505500*delt,N2);
  x15=x0+delt*(0.1*fx0-0.157178665799771*fx2+0.157178665799771*fx14);
  fx15=f(x15,theta,d+0.100000000000000000000000000000000000000000000000000000000000*delt,N2);
  x16=x0+delt*(0.181781300700095*fx0+0.675*fx1+0.34275815984719*fx2+0*fx3+0.259111214548323*fx4-0.358278966717952*fx5-1.04594895940883*fx6+0.930327845415627*fx7+1.77950959431708*fx8+0.1*fx9-0.282547569539044*fx10-0.159327350119973*fx11-0.145515894647002*fx12-0.259111214548323*fx13-0.34275815984719*fx14-0.675*fx15);
  fx16=f(x16,theta,d+delt,N2);

  x0=x0+(0.0333333333333333333333333333333333333333333333333333333333333*fx0
        +0.0250000000000000000000000000000000000000000000000000000000000*fx1
        +0.0333333333333333333333333333333333333333333333333333333333333*fx2
        +0.000000000000000000000000000000000000000000000000000000000000*fx3
        +0.0500000000000000000000000000000000000000000000000000000000000*fx4
        +0.000000000000000000000000000000000000000000000000000000000000*fx5
        +0.0400000000000000000000000000000000000000000000000000000000000*fx6
        +0.000000000000000000000000000000000000000000000000000000000000*fx7
        +0.189237478148923490158306404106012326238162346948625830327194*fx8
        +0.277429188517743176508360262560654340428504319718040836339472*fx9
        +0.277429188517743176508360262560654340428504319718040836339472*fx10
        +0.189237478148923490158306404106012326238162346948625830327194*fx11
        -0.0400000000000000000000000000000000000000000000000000000000000*fx12
        -0.0500000000000000000000000000000000000000000000000000000000000*fx13
        -0.0333333333333333333333333333333333333333333333333333333333333*fx14
        -0.0250000000000000000000000000000000000000000000000000000000000*fx15
        +0.0333333333333333333333333333333333333333333333333333333333333*fx16)*delt;
  d=d+delt;
}
'

  }
  if(!homo.res)
  {
txtB= '
  return atemp;
}

  // [[Rcpp::export]]
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N, mat delt,int N2,vec d,mat starts,int tro,int secmom)
{
  mat resss(N2,3);
  mat fx0(N2,tro);
  mat fx1(N2,tro);
  mat fx2(N2,tro);
  mat fx3(N2,tro);
  mat fx4(N2,tro);
  mat fx5(N2,tro);
  mat fx6(N2,tro);
  mat fx7(N2,tro);
  mat fx8(N2,tro);
  mat fx9(N2,tro);
  mat fx10(N2,tro);
  mat fx11(N2,tro);
  mat fx12(N2,tro);
  mat fx13(N2,tro);
  mat fx14(N2,tro);
  mat fx15(N2,tro);
  mat fx16(N2,tro);
  mat x0(N2,tro);
  mat x1(N2,tro);
  mat x2(N2,tro);
  mat x3(N2,tro);
  mat x4(N2,tro);
  mat x5(N2,tro);
  mat x6(N2,tro);
  mat x7(N2,tro);
  mat x8(N2,tro);
  mat x9(N2,tro);
  mat x10(N2,tro);
  mat x11(N2,tro);
  mat x12(N2,tro);
  mat x13(N2,tro);
  mat x14(N2,tro);
  mat x15(N2,tro);
  mat x16(N2,tro);
  x0.fill(0);
  x0.col(0)=Xs;
  x0.col(secmom-1)=Ys;
  for (int i = 1; i < N; i++)
{
  fx0=f(x0,theta,d,N2);
  x1=x0+delt%(0.1*fx0);
  fx1=f(x1,theta,d+0.100000000000000000000000000000000000000000000000000000000000*delt.col(0),N2);
  x2=x0+delt%(-0.915176561375291*fx0+1.45453440217827*fx1);
  fx2=f(x2,theta,d+0.539357840802981787532485197881302436857273449701009015505500*delt.col(0),N2);
  x3=x0+delt%(0.202259190301118*fx0+0.606777570903354*fx2);
  fx3=f(x3,theta,d+0.809036761204472681298727796821953655285910174551513523258250*delt.col(0),N2);
  x4=x0+delt%(0.184024714708644*fx0+0.197966831227192*fx2-0.0729547847313633*fx3);
  fx4=f(x4,theta,d+0.309036761204472681298727796821953655285910174551513523258250*delt.col(0),N2);
  x5=x0+delt%(0.0879007340206681*fx0+0.410459702520261*fx3+0.482713753678866*fx4);
  fx5=f(x5,theta,d+0.981074190219795268254879548310562080489056746118724882027805*delt.col(0),N2);
  x6=x0+delt%(0.085970050490246*fx0+0.330885963040722*fx3+0.48966295730945*fx4-0.0731856375070851*fx5);
  fx6=f(x6,theta,d+0.833333333333333333333333333333333333333333333333333333333333*delt.col(0),N2);
  x7=x0+delt%(0.120930449125334*fx0+0.260124675758296*fx4+0.0325402621549091*fx5-0.0595780211817361*fx6);
  fx7=f(x7,theta,d+0.354017365856802376329264185948796742115824053807373968324184*delt.col(0),N2);
  x8=x0+delt%(0.110854379580391*fx0-0.0605761488255006*fx5+0.321763705601778*fx6+0.510485725608063*fx7);
  fx8=f(x8,theta,d+0.882527661964732346425501486979669075182867844268052119663791*delt.col(0),N2);
  x9=x0+delt%(0.112054414752879*fx0-0.144942775902866*fx5-0.333269719096257*fx6+0.49926922955688*fx7+0.509504608929686*fx8);
  fx9=f(x9,theta,d+0.642615758240322548157075497020439535959501736363212695909875*delt.col(0),N2);
  x10=x0+delt%(0.113976783964186*fx0-0.0768813364203357*fx5+0.239527360324391*fx6+0.397774662368095*fx7+0.0107558956873607*fx8-0.327769124164019*fx9);
  fx10=f(x10,theta,d+0.357384241759677451842924502979560464040498263636787304090125*delt.col(0),N2);
  x11=x0+delt%(0.0798314528280196*fx0-0.0520329686800603*fx5-0.0576954146168549*fx6+0.194781915712104*fx7+0.145384923188325*fx8-0.0782942710351671*fx9-0.114503299361099*fx10);
  fx11=f(x11,theta,d+0.117472338035267653574498513020330924817132155731947880336209*delt.col(0),N2);
  x12=x0+delt%(0.985115610164857*fx0+0.330885963040722*fx3+0.48966295730945*fx4-1.37896486574844*fx5-0.861164195027636*fx6+5.78428813637537*fx7+3.28807761985104*fx8-2.38633905093136*fx9-3.25479342483644*fx10-2.16343541686423*fx11);
  fx12=f(x12,theta,d+0.833333333333333333333333333333333333333333333333333333333333*delt.col(0),N2);
  x13=x0+delt%(0.895080295771633*fx0+0.197966831227192*fx2-0.0729547847313633*fx3-0.851236239662008*fx5+0.398320112318533*fx6+3.63937263181036*fx7+1.5482287703983*fx8-2.12221714704054*fx9-1.58350398545326*fx10-1.71561608285936*fx11-0.0244036405750127*fx12);
  fx13=f(x13,theta,d+0.309036761204472681298727796821953655285910174551513523258250*delt.col(0),N2);
  x14=x0+delt%(-0.915176561375291*fx0+1.45453440217827*fx1+0*fx2+0*fx3-0.777333643644968*fx4+0*fx5-0.0910895662155176*fx6+0.0910895662155176*fx12+0.777333643644968*fx13);
  fx14=f(x14,theta,d+0.539357840802981787532485197881302436857273449701009015505500*delt.col(0),N2);
  x15=x0+delt%(0.1*fx0-0.157178665799771*fx2+0.157178665799771*fx14);
  fx15=f(x15,theta,d+0.100000000000000000000000000000000000000000000000000000000000*delt.col(0),N2);
  x16=x0+delt%(0.181781300700095*fx0+0.675*fx1+0.34275815984719*fx2+0*fx3+0.259111214548323*fx4-0.358278966717952*fx5-1.04594895940883*fx6+0.930327845415627*fx7+1.77950959431708*fx8+0.1*fx9-0.282547569539044*fx10-0.159327350119973*fx11-0.145515894647002*fx12-0.259111214548323*fx13-0.34275815984719*fx14-0.675*fx15);
  fx16=f(x16,theta,d+delt.col(0),N2);

  x0=x0+(0.0333333333333333333333333333333333333333333333333333333333333*fx0
        +0.0250000000000000000000000000000000000000000000000000000000000*fx1
        +0.0333333333333333333333333333333333333333333333333333333333333*fx2
        +0.000000000000000000000000000000000000000000000000000000000000*fx3
        +0.0500000000000000000000000000000000000000000000000000000000000*fx4
        +0.000000000000000000000000000000000000000000000000000000000000*fx5
        +0.0400000000000000000000000000000000000000000000000000000000000*fx6
        +0.000000000000000000000000000000000000000000000000000000000000*fx7
        +0.189237478148923490158306404106012326238162346948625830327194*fx8
        +0.277429188517743176508360262560654340428504319718040836339472*fx9
        +0.277429188517743176508360262560654340428504319718040836339472*fx10
        +0.189237478148923490158306404106012326238162346948625830327194*fx11
        -0.0400000000000000000000000000000000000000000000000000000000000*fx12
        -0.0500000000000000000000000000000000000000000000000000000000000*fx13
        -0.0333333333333333333333333333333333333333333333333333333333333*fx14
        -0.0250000000000000000000000000000000000000000000000000000000000*fx15
        +0.0333333333333333333333333333333333333333333333333333333333333*fx16)%delt;
 d=d+delt.col(0);
}
'
  }
}

if(state1)
{
txtC=
'
  vec det=(x0.col(1)%x0.col(3)-x0.col(4)%x0.col(4));
  resss.col(0)=-log(2*3.141592653589793)-0.5*log(abs(det))-0.5*((Xt-x0.col(0))%(Xt-x0.col(0))%x0.col(3)/det-(Xt-x0.col(0))%(Yt-x0.col(2))%x0.col(4)/det-(Xt-x0.col(0))%(Yt-x0.col(2))%x0.col(4)/det+(Yt-x0.col(2))%(Yt-x0.col(2))%x0.col(1)/det);
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
}'
}
if(state3.0)
{
txtC=
'
  vec p=(1.0/3.0) *(3*(x0.col(3)/6.0)%x0.col(1) - pow(x0.col(2)/2.0,2))/pow(x0.col(3)/6.0,2);
  vec q=(1.0/27.0)*(27*pow(x0.col(3)/6.0,2)%(x0.col(0)-Xt) - 9*(x0.col(3)/6.0)%(x0.col(2)/2.0)%x0.col(1) + 2*pow(x0.col(2)/2.0,3))/pow(x0.col(3)/6.0,3);
  vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  vec th=-(x0.col(2)/2.0)/(3*(x0.col(3)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

  vec K =x0.col(0)%th+(x0.col(1)%th%th)/2.0+(x0.col(2)%th%th%th)/6.0 +(x0.col(3)%th%th%th%th)/24.0;
  vec K1=x0.col(0)   +(x0.col(1)%th)       +(x0.col(2)%th%th)/2.0    +(x0.col(3)%th%th%th)/6.0;
  vec K2=x0.col(1)   +(x0.col(2)%th)       +(x0.col(3)%th%th)/2.0;
  vec val=-0.5*log(2*3.141592653589793*K2)+(K-th%K1);

  p=(1.0/3.0) *(3*(x0.col(7)/6.0)%x0.col(5) - pow(x0.col(6)/2.0,2))/pow(x0.col(7)/6.0,2);
  q=(1.0/27.0)*(27*pow(x0.col(7)/6.0,2)%(x0.col(4)-Yt) - 9*(x0.col(7)/6.0)%(x0.col(6)/2.0)%x0.col(5) + 2*pow(x0.col(6)/2.0,3))/pow(x0.col(7)/6.0,3);
  chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  th=-(x0.col(6)/2.0)/(3*(x0.col(7)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

  K =x0.col(4)%th+(x0.col(5)%th%th)/2.0+(x0.col(6)%th%th%th)/6.0 +(x0.col(7)%th%th%th%th)/24.0;
  K1=x0.col(4)   +(x0.col(5)%th)       +(x0.col(6)%th%th)/2.0    +(x0.col(7)%th%th%th)/6.0;
  K2=x0.col(5)   +(x0.col(6)%th)       +(x0.col(7)%th%th)/2.0;
  val=val-0.5*log(2*3.141592653589793*K2)+(K-th%K1);
  resss.col(0)=val;
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);;
}'
if(Dtype=='Normal')
{
txtC=
 '
  vec val=   -0.5*log(2*3.141592653589793*x0.col(5))-0.5*pow(Yt-x0.col(4),2)/x0.col(5);
  val=val-0.5*log(2*3.141592653589793*x0.col(1))-0.5*pow(Xt-x0.col(0),2)/x0.col(1);
  resss.col(0)=val;
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
}'
}
}
if(state3.1)
{
txtC=
'
vec p=(1.0/3.0) *(3*(x0.col(7)/6.0)%x0.col(5) - pow(x0.col(6)/2.0,2))/pow(x0.col(7)/6.0,2);
vec q=(1.0/27.0)*(27*pow(x0.col(7)/6.0,2)%(x0.col(4)-Yt) - 9*(x0.col(7)/6.0)%(x0.col(6)/2.0)%x0.col(5) + 2*pow(x0.col(6)/2.0,3))/pow(x0.col(7)/6.0,3);
vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
vec th=-(x0.col(6)/2.0)/(3*(x0.col(7)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

vec K =x0.col(4)%th+(x0.col(5)%th%th)/2.0+(x0.col(6)%th%th%th)/6.0 +(x0.col(7)%th%th%th%th)/24.0;
vec K1=x0.col(4)   +(x0.col(5)%th)       +(x0.col(6)%th%th)/2.0    +(x0.col(7)%th%th%th)/6.0;
vec K2=x0.col(5)   +(x0.col(6)%th)       +(x0.col(7)%th%th)/2.0;
vec val=-0.5*log(2*3.141592653589793*K2)+(K-th%K1);

val=val-0.5*log(2*3.141592653589793*x0.col(1))-0.5*pow(Xt-x0.col(0),2)/x0.col(1);
resss.col(0)=val;
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
  }'
if(Dtype=='Normal')
{

    txtC='
  vec val=   -0.5*log(2*3.141592653589793*x0.col(5))-0.5*pow(Yt-x0.col(4),2)/x0.col(5);
  val=val-0.5*log(2*3.141592653589793*x0.col(1))-0.5*pow(Xt-x0.col(0),2)/x0.col(1);
  resss.col(0)=val;
   List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
}'
}
}
if(state3.2)
{
txtC='
vec p=(1.0/3.0) *(3*(x0.col(3)/6.0)%x0.col(1) - pow(x0.col(2)/2.0,2))/pow(x0.col(3)/6.0,2);
vec q=(1.0/27.0)*(27*pow(x0.col(3)/6.0,2)%(x0.col(0)-Xt) - 9*(x0.col(3)/6.0)%(x0.col(2)/2.0)%x0.col(1) + 2*pow(x0.col(2)/2.0,3))/pow(x0.col(3)/6.0,3);
vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
vec th=-(x0.col(2)/2.0)/(3*(x0.col(3)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

vec K =x0.col(0)%th+(x0.col(1)%th%th)/2.0+(x0.col(2)%th%th%th)/6.0 +(x0.col(3)%th%th%th%th)/24.0;
vec K1=x0.col(0)   +(x0.col(1)%th)       +(x0.col(2)%th%th)/2.0    +(x0.col(3)%th%th%th)/6.0;
vec K2=x0.col(1)   +(x0.col(2)%th)       +(x0.col(3)%th%th)/2.0;
vec val=-0.5*log(2*3.141592653589793*K2)+(K-th%K1);

val=val-0.5*log(2*3.141592653589793*x0.col(5))-(0.5*pow(Yt-x0.col(4),2)/x0.col(5));
resss.col(0)=val;
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
}'
if(Dtype=='Normal')
{

    txtC='
  vec val=   -0.5*log(2*3.141592653589793*x0.col(5))-0.5*pow(Yt-x0.col(4),2)/x0.col(5);
  val=val-0.5*log(2*3.141592653589793*x0.col(1))-0.5*pow(Xt-x0.col(0),2)/x0.col(1);
  resss.col(0)=val;
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
}'
}
}
if(state4)
{
if(Dtype=='Saddlepoint')
{
txtC='
vec a(N2);
vec b(N2);
vec abser(N2);
abser=0.1+abser;
a.ones();
b.ones();
a=starts.col(0);
b=starts.col(1);
vec gg(N2);
vec hh(N2);
vec gg1(N2);
vec hh1(N2);
vec gg2(N2);
vec hh2(N2);
vec ar(N2);
vec br(N2);
vec anew(N2);
vec bnew(N2);
int ind=0;
while((max(abser)>0.001)&&(ind<5000))
{
gg=x0.col(0)+x0.col(1)%a+(1.0/2.0)*x0.col(2)%a%a+(1.0/6.0)*x0.col(3)%a%a%a +x0.col(8)%b +(1.0/2.0)*x0.col(9)%b%b+x0.col(10)%a%b+(1.0/6.0)*b%b%b%x0.col(12)+(1.0/2.0)*a%a%b%x0.col(13)+(1.0/2.0)*a%b%b%x0.col(11)-Xt;
hh=x0.col(4)+x0.col(5)%b+(1.0/2.0)*x0.col(6)%b%b+(1.0/6.0)*x0.col(7)%b%b%b +x0.col(8)%a +x0.col(9)%a%b+(1.0/2.0)*x0.col(10)%a%a+(1.0/2.0)*a%b%b%x0.col(12)+(1.0/6.0)*a%a%a%x0.col(13)+(1.0/2.0)*a%a%b%x0.col(11)-Yt;

gg1=x0.col(1)+x0.col(2)%a+(1.0/2.0)*x0.col(3)%a%a+x0.col(10)%b+a%b%x0.col(13)+(1.0/2.0)*b%b%x0.col(11);
gg2=x0.col(8) +x0.col(9)%b+x0.col(10)%a+(1.0/2.0)*b%b%x0.col(12)+(1.0/2.0)*a%a%x0.col(13)+a%b%x0.col(11);

hh1=x0.col(8) +x0.col(9)%b+x0.col(10)%a+(1.0/2.0)*b%b%x0.col(12)+(1.0/2.0)*a%a%x0.col(13)+a%b%x0.col(11);
hh2=x0.col(5)+x0.col(6)%b+(1.0/2.0)*x0.col(7)%b%b +x0.col(9)%a+a%b%x0.col(12)+(1.0/2.0)*a%a%x0.col(11);

anew  =a+(hh%gg2-gg%hh2)/(gg1%hh2-gg2%hh1);
bnew  =b-(hh%gg1-gg%hh1)/(gg1%hh2-gg2%hh1);
abser =(pow(anew-a,2)+pow(bnew-b,2));
a=anew;
b=bnew;
ind=ind+1;
}

resss.col(0)=-log(2*3.141592653589793)-0.5*log(gg1%hh2-gg2%hh1)+(x0.col(0)%a+x0.col(4)%b+(1.0/2.0)*x0.col(1)%a%a+(1.0/2.0)*x0.col(5)%b%b+(1.0/6.0)*x0.col(2)%a%a%a+(1.0/6.0)*x0.col(6)%b%b%b+(1.0/24.0)*x0.col(3)%a%a%a%a+(1.0/24.0)*x0.col(7)%b%b%b%b+x0.col(8)%a%b+(1.0/2.0)*x0.col(9)%a%b%b+(1.0/2.0)*x0.col(10)%a%a%b+(1.0/6.0)*a%b%b%b%x0.col(12)+(1.0/6.0)*a%a%a%b%x0.col(13)+(1.0/4.0)*a%a%b%b%x0.col(11)-a%Xt-b%Yt);
resss.col(1)=a;
resss.col(2)=b;
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
  }'
}
if(Dtype=='Edgeworth')
{
  txtC='
  vec x = (Xt-x0.col(0))/sqrt(x0.col(1));
  vec y = (Yt-x0.col(4))/sqrt(x0.col(5));

  vec rho30 =x0.col(2)/(pow(x0.col(1),1.5));
  vec rho40 =x0.col(3)/(pow(x0.col(1),2));
  vec rho03 =x0.col(6)/(pow(x0.col(5),1.5));
  vec rho04 =x0.col(7)/(pow(x0.col(5),2));
  vec rho11 =x0.col(8)/(sqrt(x0.col(1))%sqrt(x0.col(5)));
  vec rho12 =x0.col(9)/(sqrt(x0.col(1))%x0.col(5));
  vec rho21 =x0.col(10)/(x0.col(1)%sqrt(x0.col(5)));
  vec rho22 =x0.col(11)/(x0.col(1)%x0.col(5));
  vec rho13 =x0.col(12)/(sqrt(x0.col(1))%pow(x0.col(5),1.5));
  vec rho31 =x0.col(13)/(pow(x0.col(1),1.5)%sqrt(x0.col(5)));

  vec Hx1 = x;
  vec Hx2 = pow(x,2)-1;
  vec Hx3 = pow(x,3)-3*x;
  vec Hx4 = pow(x,4)-6*pow(x,2)+3;
  vec Hx5 = pow(x,5)-10*pow(x,3)+15*x;
  vec Hx6 = pow(x,6)-15*pow(x,4)+45*pow(x,2)-15;

  vec Hy1 = y;
  vec Hy2 = pow(y,2)-1;
  vec Hy3 = pow(y,3)-3*y;
  vec Hy4 = pow(y,4)-6*pow(y,2)+3;
  vec Hy5 = pow(y,5)-10*pow(y,3)+15*y;
  vec Hy6 = pow(y,6)-15*pow(y,4)+45*pow(y,2)-15;

  vec A=Hx3%rho30+3*Hx2%Hy1%rho21+3*Hx1%Hy2%rho12+Hy3%rho03;
  vec B=Hx4%rho40+4*Hx3%Hy1%rho31+6*Hx2%Hy2%rho22+4*Hx1%Hy3%rho13+Hy4%rho04;
  vec C=Hx6%pow(rho30,2)+6*Hx5%Hy1%rho21%rho30+6*Hx4%Hy2%rho12%rho30+2*Hx3%Hy3%rho03%rho30+9*Hx4%Hy2%pow(rho21,2)+18*Hx3%Hy3%rho12%rho21+6*Hx2%Hy4%rho03%rho21+9*Hx2%Hy4%pow(rho12,2)+6*Hx1%Hy5%rho03%rho12+Hy6%pow(rho03,2);

  resss.col(0)=-log(2*pi*sqrt(x0.col(1)%x0.col(5))%sqrt(1 - pow(rho11,2)))+(-(pow((Xt-x0.col(0)),2)/x0.col(1) -2*rho11%(Xt-x0.col(0))%(Yt-x0.col(4))/sqrt(x0.col(1)%x0.col(5))+pow((Yt-x0.col(4)),2)/x0.col(5))/(2*(1-pow(rho11,2))))+log(abs(1+0.1666667*A+0.04166667*B+0.01388889*C));

  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
}'
}
if(Dtype=='Normal')
{

    txtC='
   vec det=(x0.col(1)%x0.col(5)-x0.col(8)%x0.col(8));
  resss.col(0)=-log(2*3.141592653589793)-0.5*log(abs(det))-0.5*((Xt-x0.col(0))%(Xt-x0.col(0))%x0.col(5)/det-(Xt-x0.col(0))%(Yt-x0.col(4))%x0.col(8)/det-(Xt-x0.col(0))%(Yt-x0.col(4))%x0.col(8)/det+(Yt-x0.col(4))%(Yt-x0.col(4))%x0.col(1)/det);
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  return(ret);
}'
}
}

   #==============================================================================
   #                   Generate TYPE of Solution
   #==============================================================================

if(state1)
{
   # DATA RESOLUTION -------------------------------------------------------------
   if(!homo.res)
   {
     delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh)
     if(!is.null(exclude))
     {
       diffs=diff(T.seq)/mesh
       diffs[excl]=1/20/mesh
       delt=cbind(diffs,diffs,diffs,diffs,diffs)
     }
   }
   # REFERENCE MATRIX ------------------------------------------------------------
   MAT=rbind(
   c('','','(+1*a.col(0))','(+1*a.col(2))','','','','','',''),
   c('','','(+2*a.col(1))','(+2*a.col(4))','','','','','',''),
   c('','','','','(+1*a.col(0))','(+1*a.col(2))','','','',''),
   c('','','','','(+2*a.col(4))','(+2*a.col(3))','','','',''),
   c('','','(+1*a.col(4))','(+1*a.col(3))','(+1*a.col(1))','(+1*a.col(4))','','','',''))
   MAT[1,1]='(1+0*a.col(0))'        #ONES?
   MAT[3,2]='(1+0*a.col(0))'
   MAT[2,7]='(1+0*a.col(0))'
   MAT[4,10]='(1+0*a.col(0))'
   MAT[5,8]='(1+0*a.col(0))'
   MAT[5,9]='(1+0*a.col(0))'

   # HOMOGENEITY -----------------------------------------------------------------
   namess2=c('a00','b00','a10','a01','b10','b01','c00','d00','e00','f00')
   func.list2=rep(0,10)
   obs=objects(pos=1)
   for(i in 1:10)
   {
    if(sum(obs==namess2[i])){func.list2[i]=1}
   }

   func.list.timehomo=func.list2*0

   for(i in which(func.list2==1))
   {
    # which expressions vary over time
     result=eval(body(namess2[i]))
     func.list.timehomo[i]=2-(sum(diff(result)==0)==(length(result)-1))
   }
   if(any(func.list.timehomo==2)){homo=F}

   #func.list.timehomo[c(1:2,7:10)]=1    # Always set these to 1

   #BUILD ODE --------------------------------------------------------------------
   dims=rep('(',5)
   for(i in 1:5)
   {
    for(j in which(func.list2==1))
    {
      if(MAT[i,j]!='')
      {
        dims[i]=paste0(dims[i],'+(',body(namess2[j])[2],')',c('*','%')[func.list.timehomo[j]],MAT[i,j])
      }
    }
    dims[i]=paste0(dims[i],')')
   }

   if(any(dims=='()')){dims[which(dims=='()')]='(1+0*a.col(0))'}


   for(i in 1:5)
   {
      dims[i]=paste0(paste0(paste0('   atemp.col(',i-1,')='),dims[i]),';')
   }

   # WRIGHT AND SOURCE -----------------------------------------------------------
   txt.full=paste(txtA,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],txtB,txtC)
   type.sol ="                  GENERALIZED LINEAR DIFFUSON"
}
if(any(c(state3.0,state3.1,state3.2)))
{
   # DATA RESOLUTION -------------------------------------------------------------
   if(!homo.res)
   {
     delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh)
     if(!is.null(exclude))
     {
       diffs=diff(T.seq)/mesh
       diffs[excl]=1/20/mesh
       delt=cbind(diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs)
     }
   }
   # REFERENCE MATRIX ------------------------------------------------------------
   MAT=rbind(
    c("()","(+1*a.col(0))","(+1*a.col(0)%a.col(0)+1*a.col(1))","(+1*a.col(4))","(+1*a.col(4)%a.col(4)+1*a.col(5))","(+1*a.col(8)+1*a.col(0)%a.col(4))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","(+2*a.col(1))","(+2*a.col(0)%a.col(1)+2*a.col(1)%a.col(0)+2*a.col(2))","(+2*a.col(8))","(+2*a.col(4)%a.col(8)+2*a.col(8)%a.col(4)+2*a.col(9))","(+2*a.col(10)+2*a.col(0)%a.col(8)+2*a.col(1)%a.col(4))","()","()","()","()","()","()","()","(+1*a.col(0))","(+1*a.col(1)+1*a.col(0)%a.col(0))","(+1*a.col(4))","(+1*a.col(5)+1*a.col(4)%a.col(4))","(+1*a.col(8)+1*a.col(4)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","(+3*a.col(2))","(+3*a.col(0)%a.col(2)+6*a.col(1)%a.col(1)+3*a.col(2)%a.col(0)+3*a.col(3))","(+3*a.col(10))","(+3*a.col(4)%a.col(10)+6*a.col(8)%a.col(8)+3*a.col(10)%a.col(4)+3*a.col(11))","(+3*a.col(13)+3*a.col(0)%a.col(10)+6*a.col(1)%a.col(8)+3*a.col(2)%a.col(4))","()","()","()","()","()","()","()","(+3*a.col(1))","(+3*a.col(2)+3*a.col(0)%a.col(1)+3*a.col(1)%a.col(0))","(+3*a.col(8))","(+3*a.col(9)+3*a.col(4)%a.col(8)+3*a.col(8)%a.col(4))","(+3*a.col(10)+3*a.col(4)%a.col(1)+3*a.col(8)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","(+4*a.col(3))","(+4*a.col(0)%a.col(3)+12*a.col(1)%a.col(2)+12*a.col(2)%a.col(1)+4*a.col(3)%a.col(0))","(+4*a.col(13))","(+4*a.col(4)%a.col(13)+12*a.col(8)%a.col(10)+12*a.col(10)%a.col(8)+4*a.col(13)%a.col(4))","(+4*a.col(0)%a.col(13)+12*a.col(1)%a.col(10)+12*a.col(2)%a.col(8)+4*a.col(3)%a.col(4))","()","()","()","()","()","()","()","(+6*a.col(2))","(+6*a.col(3)+6*a.col(0)%a.col(2)+12*a.col(1)%a.col(1)+6*a.col(2)%a.col(0))","(+6*a.col(10))","(+6*a.col(11)+6*a.col(4)%a.col(10)+12*a.col(8)%a.col(8)+6*a.col(10)%a.col(4))","(+6*a.col(13)+6*a.col(4)%a.col(2)+12*a.col(8)%a.col(1)+6*a.col(10)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","()","()","()","()","()","()","(+1*a.col(0))","(+1*a.col(0)%a.col(0)+1*a.col(1))","(+1*a.col(4))","(+1*a.col(4)%a.col(4)+1*a.col(5))","(+1*a.col(8)+1*a.col(4)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","()","()","()","()","()","()","(+2*a.col(8))","(+2*a.col(0)%a.col(8)+2*a.col(8)%a.col(0)+2*a.col(10))","(+2*a.col(5))","(+2*a.col(4)%a.col(5)+2*a.col(5)%a.col(4)+2*a.col(6))","(+2*a.col(9)+2*a.col(4)%a.col(8)+2*a.col(5)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","(+1*a.col(0))","(+1*a.col(1)+1*a.col(0)%a.col(0))","(+1*a.col(4))","(+1*a.col(5)+1*a.col(4)%a.col(4))","(+1*a.col(8)+1*a.col(4)%a.col(0))"
    ),c("()","()","()","()","()","()","()","(+3*a.col(9))","(+3*a.col(0)%a.col(9)+6*a.col(8)%a.col(8)+3*a.col(9)%a.col(0)+3*a.col(11))","(+3*a.col(6))","(+3*a.col(4)%a.col(6)+6*a.col(5)%a.col(5)+3*a.col(6)%a.col(4)+3*a.col(7))","(+3*a.col(12)+3*a.col(4)%a.col(9)+6*a.col(5)%a.col(8)+3*a.col(6)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","(+3*a.col(8))","(+3*a.col(10)+3*a.col(0)%a.col(8)+3*a.col(8)%a.col(0))","(+3*a.col(5))","(+3*a.col(6)+3*a.col(4)%a.col(5)+3*a.col(5)%a.col(4))","(+3*a.col(9)+3*a.col(4)%a.col(8)+3*a.col(5)%a.col(0))"
    ),c("()","()","()","()","()","()","()","(+4*a.col(12))","(+4*a.col(0)%a.col(12)+12*a.col(8)%a.col(9)+12*a.col(9)%a.col(8)+4*a.col(12)%a.col(0))","(+4*a.col(7))","(+4*a.col(4)%a.col(7)+12*a.col(5)%a.col(6)+12*a.col(6)%a.col(5)+4*a.col(7)%a.col(4))","(+4*a.col(4)%a.col(12)+12*a.col(5)%a.col(9)+12*a.col(6)%a.col(8)+4*a.col(7)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","(+6*a.col(9))","(+6*a.col(11)+6*a.col(0)%a.col(9)+12*a.col(8)%a.col(8)+6*a.col(9)%a.col(0))","(+6*a.col(6))","(+6*a.col(7)+6*a.col(4)%a.col(6)+12*a.col(5)%a.col(5)+6*a.col(6)%a.col(4))","(+6*a.col(12)+6*a.col(4)%a.col(9)+12*a.col(5)%a.col(8)+6*a.col(6)%a.col(0))"
    ),c("()","(+1*a.col(8))","(+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0)+1*a.col(10))","(+1*a.col(5))","(+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4)+1*a.col(6))","(+1*a.col(9)+1*a.col(0)%a.col(5)+1*a.col(8)%a.col(4))","()","(+1*a.col(1))","(+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0)+1*a.col(2))","(+1*a.col(8))","(+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4)+1*a.col(9))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))","()","()","()","()","()","()","()","(+0.5*a.col(0))","(+0.5*a.col(1)+0.5*a.col(0)%a.col(0))","(+0.5*a.col(4))","(+0.5*a.col(5)+0.5*a.col(4)%a.col(4))","(+0.5*a.col(8)+0.5*a.col(4)%a.col(0))","()","(+0.5*a.col(0))","(+0.5*a.col(1)+0.5*a.col(0)%a.col(0))","(+0.5*a.col(4))","(+0.5*a.col(5)+0.5*a.col(4)%a.col(4))","(+0.5*a.col(8)+0.5*a.col(4)%a.col(0))","()","()","()","()","()","()"
    ),c("()","(+1*a.col(9))","(+1*a.col(0)%a.col(9)+2*a.col(8)%a.col(8)+1*a.col(9)%a.col(0)+1*a.col(11))","(+1*a.col(6))","(+1*a.col(4)%a.col(6)+2*a.col(5)%a.col(5)+1*a.col(6)%a.col(4)+1*a.col(7))","(+1*a.col(12)+1*a.col(0)%a.col(6)+2*a.col(8)%a.col(5)+1*a.col(9)%a.col(4))","()","(+2*a.col(10))","(+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0)+2*a.col(13))","(+2*a.col(9))","(+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4)+2*a.col(12))","(+2*a.col(11)+2*a.col(4)%a.col(10)+2*a.col(5)%a.col(1)+2*a.col(8)%a.col(8)+2*a.col(9)%a.col(0))","()","()","()","()","()","()","()","(+1*a.col(8))","(+1*a.col(10)+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0))","(+1*a.col(5))","(+1*a.col(6)+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(5)%a.col(0))","()","(+1*a.col(8))","(+1*a.col(10)+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0))","(+1*a.col(5))","(+1*a.col(6)+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(5)%a.col(0))","()","(+1*a.col(1))","(+1*a.col(2)+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0))","(+1*a.col(8))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))"
    ),c("()","(+2*a.col(10))","(+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0)+2*a.col(13))","(+2*a.col(9))","(+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4)+2*a.col(12))","(+2*a.col(11)+2*a.col(0)%a.col(9)+2*a.col(8)%a.col(8)+2*a.col(1)%a.col(5)+2*a.col(10)%a.col(4))","()","(+1*a.col(2))","(+1*a.col(0)%a.col(2)+2*a.col(1)%a.col(1)+1*a.col(2)%a.col(0)+1*a.col(3))","(+1*a.col(10))","(+1*a.col(4)%a.col(10)+2*a.col(8)%a.col(8)+1*a.col(10)%a.col(4)+1*a.col(11))","(+1*a.col(13)+1*a.col(4)%a.col(2)+2*a.col(8)%a.col(1)+1*a.col(10)%a.col(0))","()","(+1*a.col(8))","(+1*a.col(10)+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0))","(+1*a.col(5))","(+1*a.col(6)+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(5)%a.col(0))","()","(+1*a.col(1))","(+1*a.col(2)+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0))","(+1*a.col(8))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))","()","(+1*a.col(1))","(+1*a.col(2)+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0))","(+1*a.col(8))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))","()","()","()","()","()","()"
    ),c("()","(+2*a.col(11))","(+2*a.col(0)%a.col(11)+4*a.col(8)%a.col(10)+2*a.col(9)%a.col(1)+2*a.col(1)%a.col(9)+4*a.col(10)%a.col(8)+2*a.col(11)%a.col(0))","(+2*a.col(12))","(+2*a.col(4)%a.col(12)+4*a.col(5)%a.col(9)+2*a.col(6)%a.col(8)+2*a.col(8)%a.col(6)+4*a.col(9)%a.col(5)+2*a.col(12)%a.col(4))","(+2*a.col(0)%a.col(12)+4*a.col(8)%a.col(9)+2*a.col(9)%a.col(8)+2*a.col(1)%a.col(6)+4*a.col(10)%a.col(5)+2*a.col(11)%a.col(4))","()","(+2*a.col(13))","(+2*a.col(0)%a.col(13)+2*a.col(8)%a.col(2)+4*a.col(1)%a.col(10)+4*a.col(10)%a.col(1)+2*a.col(2)%a.col(8)+2*a.col(13)%a.col(0))","(+2*a.col(11))","(+2*a.col(4)%a.col(11)+2*a.col(5)%a.col(10)+4*a.col(8)%a.col(9)+4*a.col(9)%a.col(8)+2*a.col(10)%a.col(5)+2*a.col(11)%a.col(4))","(+2*a.col(4)%a.col(13)+2*a.col(5)%a.col(2)+4*a.col(8)%a.col(10)+4*a.col(9)%a.col(1)+2*a.col(10)%a.col(8)+2*a.col(11)%a.col(0))","()","(+1*a.col(9))","(+1*a.col(11)+1*a.col(0)%a.col(9)+2*a.col(8)%a.col(8)+1*a.col(9)%a.col(0))","(+1*a.col(6))","(+1*a.col(7)+1*a.col(4)%a.col(6)+2*a.col(5)%a.col(5)+1*a.col(6)%a.col(4))","(+1*a.col(12)+1*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+1*a.col(6)%a.col(0))","()","(+2*a.col(10))","(+2*a.col(13)+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0))","(+2*a.col(9))","(+2*a.col(12)+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4))","(+2*a.col(11)+2*a.col(4)%a.col(10)+2*a.col(5)%a.col(1)+2*a.col(8)%a.col(8)+2*a.col(9)%a.col(0))","()","(+2*a.col(10))","(+2*a.col(13)+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0))","(+2*a.col(9))","(+2*a.col(12)+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4))","(+2*a.col(11)+2*a.col(4)%a.col(10)+2*a.col(5)%a.col(1)+2*a.col(8)%a.col(8)+2*a.col(9)%a.col(0))","()","(+1*a.col(2))","(+1*a.col(3)+1*a.col(0)%a.col(2)+2*a.col(1)%a.col(1)+1*a.col(2)%a.col(0))","(+1*a.col(10))","(+1*a.col(11)+1*a.col(4)%a.col(10)+2*a.col(8)%a.col(8)+1*a.col(10)%a.col(4))","(+1*a.col(13)+1*a.col(4)%a.col(2)+2*a.col(8)%a.col(1)+1*a.col(10)%a.col(0))"
    ),c("()","(+1*a.col(12))","(+1*a.col(0)%a.col(12)+3*a.col(8)%a.col(9)+3*a.col(9)%a.col(8)+1*a.col(12)%a.col(0))","(+1*a.col(7))","(+1*a.col(4)%a.col(7)+3*a.col(5)%a.col(6)+3*a.col(6)%a.col(5)+1*a.col(7)%a.col(4))","(+1*a.col(0)%a.col(7)+3*a.col(8)%a.col(6)+3*a.col(9)%a.col(5)+1*a.col(12)%a.col(4))","()","(+3*a.col(11))","(+3*a.col(0)%a.col(11)+6*a.col(8)%a.col(10)+3*a.col(9)%a.col(1)+3*a.col(1)%a.col(9)+6*a.col(10)%a.col(8)+3*a.col(11)%a.col(0))","(+3*a.col(12))","(+3*a.col(4)%a.col(12)+6*a.col(5)%a.col(9)+3*a.col(6)%a.col(8)+3*a.col(8)%a.col(6)+6*a.col(9)%a.col(5)+3*a.col(12)%a.col(4))","(+3*a.col(4)%a.col(11)+6*a.col(5)%a.col(10)+3*a.col(6)%a.col(1)+3*a.col(8)%a.col(9)+6*a.col(9)%a.col(8)+3*a.col(12)%a.col(0))","()","()","()","()","()","()","()","(+1.5*a.col(9))","(+1.5*a.col(11)+1.5*a.col(0)%a.col(9)+3*a.col(8)%a.col(8)+1.5*a.col(9)%a.col(0))","(+1.5*a.col(6))","(+1.5*a.col(7)+1.5*a.col(4)%a.col(6)+3*a.col(5)%a.col(5)+1.5*a.col(6)%a.col(4))","(+1.5*a.col(12)+1.5*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+1.5*a.col(6)%a.col(0))","()","(+1.5*a.col(9))","(+1.5*a.col(11)+1.5*a.col(0)%a.col(9)+3*a.col(8)%a.col(8)+1.5*a.col(9)%a.col(0))","(+1.5*a.col(6))","(+1.5*a.col(7)+1.5*a.col(4)%a.col(6)+3*a.col(5)%a.col(5)+1.5*a.col(6)%a.col(4))","(+1.5*a.col(12)+1.5*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+1.5*a.col(6)%a.col(0))","()","(+3*a.col(10))","(+3*a.col(13)+3*a.col(0)%a.col(10)+3*a.col(8)%a.col(1)+3*a.col(1)%a.col(8)+3*a.col(10)%a.col(0))","(+3*a.col(9))","(+3*a.col(12)+3*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+3*a.col(8)%a.col(5)+3*a.col(9)%a.col(4))","(+3*a.col(11)+3*a.col(4)%a.col(10)+3*a.col(5)%a.col(1)+3*a.col(8)%a.col(8)+3*a.col(9)%a.col(0))"
    ),c("()","(+3*a.col(13))","(+3*a.col(0)%a.col(13)+3*a.col(8)%a.col(2)+6*a.col(1)%a.col(10)+6*a.col(10)%a.col(1)+3*a.col(2)%a.col(8)+3*a.col(13)%a.col(0))","(+3*a.col(11))","(+3*a.col(4)%a.col(11)+3*a.col(5)%a.col(10)+6*a.col(8)%a.col(9)+6*a.col(9)%a.col(8)+3*a.col(10)%a.col(5)+3*a.col(11)%a.col(4))","(+3*a.col(0)%a.col(11)+3*a.col(8)%a.col(10)+6*a.col(1)%a.col(9)+6*a.col(10)%a.col(8)+3*a.col(2)%a.col(5)+3*a.col(13)%a.col(4))","()","(+1*a.col(3))","(+1*a.col(0)%a.col(3)+3*a.col(1)%a.col(2)+3*a.col(2)%a.col(1)+1*a.col(3)%a.col(0))","(+1*a.col(13))","(+1*a.col(4)%a.col(13)+3*a.col(8)%a.col(10)+3*a.col(10)%a.col(8)+1*a.col(13)%a.col(4))","(+1*a.col(4)%a.col(3)+3*a.col(8)%a.col(2)+3*a.col(10)%a.col(1)+1*a.col(13)%a.col(0))","()","(+3*a.col(10))","(+3*a.col(13)+3*a.col(0)%a.col(10)+3*a.col(8)%a.col(1)+3*a.col(1)%a.col(8)+3*a.col(10)%a.col(0))","(+3*a.col(9))","(+3*a.col(12)+3*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+3*a.col(8)%a.col(5)+3*a.col(9)%a.col(4))","(+3*a.col(11)+3*a.col(4)%a.col(10)+3*a.col(5)%a.col(1)+3*a.col(8)%a.col(8)+3*a.col(9)%a.col(0))","()","(+1.5*a.col(2))","(+1.5*a.col(3)+1.5*a.col(0)%a.col(2)+3*a.col(1)%a.col(1)+1.5*a.col(2)%a.col(0))","(+1.5*a.col(10))","(+1.5*a.col(11)+1.5*a.col(4)%a.col(10)+3*a.col(8)%a.col(8)+1.5*a.col(10)%a.col(4))","(+1.5*a.col(13)+1.5*a.col(4)%a.col(2)+3*a.col(8)%a.col(1)+1.5*a.col(10)%a.col(0))","()","(+1.5*a.col(2))","(+1.5*a.col(3)+1.5*a.col(0)%a.col(2)+3*a.col(1)%a.col(1)+1.5*a.col(2)%a.col(0))","(+1.5*a.col(10))","(+1.5*a.col(11)+1.5*a.col(4)%a.col(10)+3*a.col(8)%a.col(8)+1.5*a.col(10)%a.col(4))","(+1.5*a.col(13)+1.5*a.col(4)%a.col(2)+3*a.col(8)%a.col(1)+1.5*a.col(10)%a.col(0))","()","()","()","()","()","()"
    ))

    MAT[1,1]='(1+0*a.col(0))'           # ONES(N2)
    MAT[5,7]='(1+0*a.col(0))'
    MAT[2,13]='(1+0*a.col(0))'
    MAT[9,19]='(1+0*a.col(0))'
    MAT[9,25]='(1+0*a.col(0))'
    MAT[6,31]='(1+0*a.col(0))'

   # HOMOGENEITY -----------------------------------------------------------------
   func.list.timehomo=func.list*0
   for(i in which(func.list==1))
   {
    # which expressions vary over time
     result=eval(body(namess[i]))
     func.list.timehomo[i]=2-(sum(diff(result)==0)==(length(result)-1))
   }
   if(any(func.list.timehomo==2)){homo=F}
   #func.list.timehomo[c(1+(0:5)*6)]=1    # Always set these to 1


   # BUILD ODE -------------------------------------------------------------------
   dims=rep('(',8)
   for(i in 1:8)
   {
     for(j in which(func.list==1))
     {
         if(MAT[i,j]!="()")
         {
           dims[i]=paste0(dims[i],'+(',body(namess[j])[2],')',c('*','%')[func.list.timehomo[j]],MAT[i,j])
         }
     }
     dims[i]=paste0(dims[i],')')
   }
    if(any(dims=='()')){dims[which(dims=='()')]='(1+0*a.col(0))'}
   for(i in 1:8)
   {
        dims[i]=paste0(paste0(paste0('   atemp.col(',i-1,')='),dims[i]),';')
   }

   # WRIGHT AND SOURCE -----------------------------------------------------------
   txt.full=paste(txtA,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],'\n',dims[6],'\n',dims[7],'\n',dims[8],txtB,txtC)
   type.sol ="                  GENERALIZED QUADRATIC DIFFUSON"

}

if(state4)
{
   # DATA RESOLUTION -------------------------------------------------------------
   if(!homo.res)
   {
     delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,
                diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh)
     if(!is.null(exclude))
     {
       diffs=diff(T.seq)/mesh
       diffs[excl]=1/20/mesh
       delt=cbind(diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs,diffs)
     }
   }
   # REFERENCE MATRIX ------------------------------------------------------------
   MAT=rbind(
    c("()","(+1*a.col(0))","(+1*a.col(0)%a.col(0)+1*a.col(1))","(+1*a.col(4))","(+1*a.col(4)%a.col(4)+1*a.col(5))","(+1*a.col(8)+1*a.col(0)%a.col(4))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","(+2*a.col(1))","(+2*a.col(0)%a.col(1)+2*a.col(1)%a.col(0)+2*a.col(2))","(+2*a.col(8))","(+2*a.col(4)%a.col(8)+2*a.col(8)%a.col(4)+2*a.col(9))","(+2*a.col(10)+2*a.col(0)%a.col(8)+2*a.col(1)%a.col(4))","()","()","()","()","()","()","()","(+1*a.col(0))","(+1*a.col(1)+1*a.col(0)%a.col(0))","(+1*a.col(4))","(+1*a.col(5)+1*a.col(4)%a.col(4))","(+1*a.col(8)+1*a.col(4)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","(+3*a.col(2))","(+3*a.col(0)%a.col(2)+6*a.col(1)%a.col(1)+3*a.col(2)%a.col(0)+3*a.col(3))","(+3*a.col(10))","(+3*a.col(4)%a.col(10)+6*a.col(8)%a.col(8)+3*a.col(10)%a.col(4)+3*a.col(11))","(+3*a.col(13)+3*a.col(0)%a.col(10)+6*a.col(1)%a.col(8)+3*a.col(2)%a.col(4))","()","()","()","()","()","()","()","(+3*a.col(1))","(+3*a.col(2)+3*a.col(0)%a.col(1)+3*a.col(1)%a.col(0))","(+3*a.col(8))","(+3*a.col(9)+3*a.col(4)%a.col(8)+3*a.col(8)%a.col(4))","(+3*a.col(10)+3*a.col(4)%a.col(1)+3*a.col(8)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","(+4*a.col(3))","(+4*a.col(0)%a.col(3)+12*a.col(1)%a.col(2)+12*a.col(2)%a.col(1)+4*a.col(3)%a.col(0))","(+4*a.col(13))","(+4*a.col(4)%a.col(13)+12*a.col(8)%a.col(10)+12*a.col(10)%a.col(8)+4*a.col(13)%a.col(4))","(+4*a.col(0)%a.col(13)+12*a.col(1)%a.col(10)+12*a.col(2)%a.col(8)+4*a.col(3)%a.col(4))","()","()","()","()","()","()","()","(+6*a.col(2))","(+6*a.col(3)+6*a.col(0)%a.col(2)+12*a.col(1)%a.col(1)+6*a.col(2)%a.col(0))","(+6*a.col(10))","(+6*a.col(11)+6*a.col(4)%a.col(10)+12*a.col(8)%a.col(8)+6*a.col(10)%a.col(4))","(+6*a.col(13)+6*a.col(4)%a.col(2)+12*a.col(8)%a.col(1)+6*a.col(10)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","()","()","()","()","()","()","(+1*a.col(0))","(+1*a.col(0)%a.col(0)+1*a.col(1))","(+1*a.col(4))","(+1*a.col(4)%a.col(4)+1*a.col(5))","(+1*a.col(8)+1*a.col(4)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()"
    ),c("()","()","()","()","()","()","()","(+2*a.col(8))","(+2*a.col(0)%a.col(8)+2*a.col(8)%a.col(0)+2*a.col(10))","(+2*a.col(5))","(+2*a.col(4)%a.col(5)+2*a.col(5)%a.col(4)+2*a.col(6))","(+2*a.col(9)+2*a.col(4)%a.col(8)+2*a.col(5)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","(+1*a.col(0))","(+1*a.col(1)+1*a.col(0)%a.col(0))","(+1*a.col(4))","(+1*a.col(5)+1*a.col(4)%a.col(4))","(+1*a.col(8)+1*a.col(4)%a.col(0))"
    ),c("()","()","()","()","()","()","()","(+3*a.col(9))","(+3*a.col(0)%a.col(9)+6*a.col(8)%a.col(8)+3*a.col(9)%a.col(0)+3*a.col(11))","(+3*a.col(6))","(+3*a.col(4)%a.col(6)+6*a.col(5)%a.col(5)+3*a.col(6)%a.col(4)+3*a.col(7))","(+3*a.col(12)+3*a.col(4)%a.col(9)+6*a.col(5)%a.col(8)+3*a.col(6)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","(+3*a.col(8))","(+3*a.col(10)+3*a.col(0)%a.col(8)+3*a.col(8)%a.col(0))","(+3*a.col(5))","(+3*a.col(6)+3*a.col(4)%a.col(5)+3*a.col(5)%a.col(4))","(+3*a.col(9)+3*a.col(4)%a.col(8)+3*a.col(5)%a.col(0))"
    ),c("()","()","()","()","()","()","()","(+4*a.col(12))","(+4*a.col(0)%a.col(12)+12*a.col(8)%a.col(9)+12*a.col(9)%a.col(8)+4*a.col(12)%a.col(0))","(+4*a.col(7))","(+4*a.col(4)%a.col(7)+12*a.col(5)%a.col(6)+12*a.col(6)%a.col(5)+4*a.col(7)%a.col(4))","(+4*a.col(4)%a.col(12)+12*a.col(5)%a.col(9)+12*a.col(6)%a.col(8)+4*a.col(7)%a.col(0))","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","()","(+6*a.col(9))","(+6*a.col(11)+6*a.col(0)%a.col(9)+12*a.col(8)%a.col(8)+6*a.col(9)%a.col(0))","(+6*a.col(6))","(+6*a.col(7)+6*a.col(4)%a.col(6)+12*a.col(5)%a.col(5)+6*a.col(6)%a.col(4))","(+6*a.col(12)+6*a.col(4)%a.col(9)+12*a.col(5)%a.col(8)+6*a.col(6)%a.col(0))"
    ),c("()","(+1*a.col(8))","(+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0)+1*a.col(10))","(+1*a.col(5))","(+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4)+1*a.col(6))","(+1*a.col(9)+1*a.col(0)%a.col(5)+1*a.col(8)%a.col(4))","()","(+1*a.col(1))","(+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0)+1*a.col(2))","(+1*a.col(8))","(+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4)+1*a.col(9))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))","()","()","()","()","()","()","()","(+0.5*a.col(0))","(+0.5*a.col(1)+0.5*a.col(0)%a.col(0))","(+0.5*a.col(4))","(+0.5*a.col(5)+0.5*a.col(4)%a.col(4))","(+0.5*a.col(8)+0.5*a.col(4)%a.col(0))","()","(+0.5*a.col(0))","(+0.5*a.col(1)+0.5*a.col(0)%a.col(0))","(+0.5*a.col(4))","(+0.5*a.col(5)+0.5*a.col(4)%a.col(4))","(+0.5*a.col(8)+0.5*a.col(4)%a.col(0))","()","()","()","()","()","()"
    ),c("()","(+1*a.col(9))","(+1*a.col(0)%a.col(9)+2*a.col(8)%a.col(8)+1*a.col(9)%a.col(0)+1*a.col(11))","(+1*a.col(6))","(+1*a.col(4)%a.col(6)+2*a.col(5)%a.col(5)+1*a.col(6)%a.col(4)+1*a.col(7))","(+1*a.col(12)+1*a.col(0)%a.col(6)+2*a.col(8)%a.col(5)+1*a.col(9)%a.col(4))","()","(+2*a.col(10))","(+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0)+2*a.col(13))","(+2*a.col(9))","(+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4)+2*a.col(12))","(+2*a.col(11)+2*a.col(4)%a.col(10)+2*a.col(5)%a.col(1)+2*a.col(8)%a.col(8)+2*a.col(9)%a.col(0))","()","()","()","()","()","()","()","(+1*a.col(8))","(+1*a.col(10)+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0))","(+1*a.col(5))","(+1*a.col(6)+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(5)%a.col(0))","()","(+1*a.col(8))","(+1*a.col(10)+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0))","(+1*a.col(5))","(+1*a.col(6)+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(5)%a.col(0))","()","(+1*a.col(1))","(+1*a.col(2)+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0))","(+1*a.col(8))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))"
    ),c("()","(+2*a.col(10))","(+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0)+2*a.col(13))","(+2*a.col(9))","(+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4)+2*a.col(12))","(+2*a.col(11)+2*a.col(0)%a.col(9)+2*a.col(8)%a.col(8)+2*a.col(1)%a.col(5)+2*a.col(10)%a.col(4))","()","(+1*a.col(2))","(+1*a.col(0)%a.col(2)+2*a.col(1)%a.col(1)+1*a.col(2)%a.col(0)+1*a.col(3))","(+1*a.col(10))","(+1*a.col(4)%a.col(10)+2*a.col(8)%a.col(8)+1*a.col(10)%a.col(4)+1*a.col(11))","(+1*a.col(13)+1*a.col(4)%a.col(2)+2*a.col(8)%a.col(1)+1*a.col(10)%a.col(0))","()","(+1*a.col(8))","(+1*a.col(10)+1*a.col(0)%a.col(8)+1*a.col(8)%a.col(0))","(+1*a.col(5))","(+1*a.col(6)+1*a.col(4)%a.col(5)+1*a.col(5)%a.col(4))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(5)%a.col(0))","()","(+1*a.col(1))","(+1*a.col(2)+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0))","(+1*a.col(8))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))","()","(+1*a.col(1))","(+1*a.col(2)+1*a.col(0)%a.col(1)+1*a.col(1)%a.col(0))","(+1*a.col(8))","(+1*a.col(9)+1*a.col(4)%a.col(8)+1*a.col(8)%a.col(4))","(+1*a.col(10)+1*a.col(4)%a.col(1)+1*a.col(8)%a.col(0))","()","()","()","()","()","()"
    ),c("()","(+2*a.col(11))","(+2*a.col(0)%a.col(11)+4*a.col(8)%a.col(10)+2*a.col(9)%a.col(1)+2*a.col(1)%a.col(9)+4*a.col(10)%a.col(8)+2*a.col(11)%a.col(0))","(+2*a.col(12))","(+2*a.col(4)%a.col(12)+4*a.col(5)%a.col(9)+2*a.col(6)%a.col(8)+2*a.col(8)%a.col(6)+4*a.col(9)%a.col(5)+2*a.col(12)%a.col(4))","(+2*a.col(0)%a.col(12)+4*a.col(8)%a.col(9)+2*a.col(9)%a.col(8)+2*a.col(1)%a.col(6)+4*a.col(10)%a.col(5)+2*a.col(11)%a.col(4))","()","(+2*a.col(13))","(+2*a.col(0)%a.col(13)+2*a.col(8)%a.col(2)+4*a.col(1)%a.col(10)+4*a.col(10)%a.col(1)+2*a.col(2)%a.col(8)+2*a.col(13)%a.col(0))","(+2*a.col(11))","(+2*a.col(4)%a.col(11)+2*a.col(5)%a.col(10)+4*a.col(8)%a.col(9)+4*a.col(9)%a.col(8)+2*a.col(10)%a.col(5)+2*a.col(11)%a.col(4))","(+2*a.col(4)%a.col(13)+2*a.col(5)%a.col(2)+4*a.col(8)%a.col(10)+4*a.col(9)%a.col(1)+2*a.col(10)%a.col(8)+2*a.col(11)%a.col(0))","()","(+1*a.col(9))","(+1*a.col(11)+1*a.col(0)%a.col(9)+2*a.col(8)%a.col(8)+1*a.col(9)%a.col(0))","(+1*a.col(6))","(+1*a.col(7)+1*a.col(4)%a.col(6)+2*a.col(5)%a.col(5)+1*a.col(6)%a.col(4))","(+1*a.col(12)+1*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+1*a.col(6)%a.col(0))","()","(+2*a.col(10))","(+2*a.col(13)+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0))","(+2*a.col(9))","(+2*a.col(12)+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4))","(+2*a.col(11)+2*a.col(4)%a.col(10)+2*a.col(5)%a.col(1)+2*a.col(8)%a.col(8)+2*a.col(9)%a.col(0))","()","(+2*a.col(10))","(+2*a.col(13)+2*a.col(0)%a.col(10)+2*a.col(8)%a.col(1)+2*a.col(1)%a.col(8)+2*a.col(10)%a.col(0))","(+2*a.col(9))","(+2*a.col(12)+2*a.col(4)%a.col(9)+2*a.col(5)%a.col(8)+2*a.col(8)%a.col(5)+2*a.col(9)%a.col(4))","(+2*a.col(11)+2*a.col(4)%a.col(10)+2*a.col(5)%a.col(1)+2*a.col(8)%a.col(8)+2*a.col(9)%a.col(0))","()","(+1*a.col(2))","(+1*a.col(3)+1*a.col(0)%a.col(2)+2*a.col(1)%a.col(1)+1*a.col(2)%a.col(0))","(+1*a.col(10))","(+1*a.col(11)+1*a.col(4)%a.col(10)+2*a.col(8)%a.col(8)+1*a.col(10)%a.col(4))","(+1*a.col(13)+1*a.col(4)%a.col(2)+2*a.col(8)%a.col(1)+1*a.col(10)%a.col(0))"
    ),c("()","(+1*a.col(12))","(+1*a.col(0)%a.col(12)+3*a.col(8)%a.col(9)+3*a.col(9)%a.col(8)+1*a.col(12)%a.col(0))","(+1*a.col(7))","(+1*a.col(4)%a.col(7)+3*a.col(5)%a.col(6)+3*a.col(6)%a.col(5)+1*a.col(7)%a.col(4))","(+1*a.col(0)%a.col(7)+3*a.col(8)%a.col(6)+3*a.col(9)%a.col(5)+1*a.col(12)%a.col(4))","()","(+3*a.col(11))","(+3*a.col(0)%a.col(11)+6*a.col(8)%a.col(10)+3*a.col(9)%a.col(1)+3*a.col(1)%a.col(9)+6*a.col(10)%a.col(8)+3*a.col(11)%a.col(0))","(+3*a.col(12))","(+3*a.col(4)%a.col(12)+6*a.col(5)%a.col(9)+3*a.col(6)%a.col(8)+3*a.col(8)%a.col(6)+6*a.col(9)%a.col(5)+3*a.col(12)%a.col(4))","(+3*a.col(4)%a.col(11)+6*a.col(5)%a.col(10)+3*a.col(6)%a.col(1)+3*a.col(8)%a.col(9)+6*a.col(9)%a.col(8)+3*a.col(12)%a.col(0))","()","()","()","()","()","()","()","(+1.5*a.col(9))","(+1.5*a.col(11)+1.5*a.col(0)%a.col(9)+3*a.col(8)%a.col(8)+1.5*a.col(9)%a.col(0))","(+1.5*a.col(6))","(+1.5*a.col(7)+1.5*a.col(4)%a.col(6)+3*a.col(5)%a.col(5)+1.5*a.col(6)%a.col(4))","(+1.5*a.col(12)+1.5*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+1.5*a.col(6)%a.col(0))","()","(+1.5*a.col(9))","(+1.5*a.col(11)+1.5*a.col(0)%a.col(9)+3*a.col(8)%a.col(8)+1.5*a.col(9)%a.col(0))","(+1.5*a.col(6))","(+1.5*a.col(7)+1.5*a.col(4)%a.col(6)+3*a.col(5)%a.col(5)+1.5*a.col(6)%a.col(4))","(+1.5*a.col(12)+1.5*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+1.5*a.col(6)%a.col(0))","()","(+3*a.col(10))","(+3*a.col(13)+3*a.col(0)%a.col(10)+3*a.col(8)%a.col(1)+3*a.col(1)%a.col(8)+3*a.col(10)%a.col(0))","(+3*a.col(9))","(+3*a.col(12)+3*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+3*a.col(8)%a.col(5)+3*a.col(9)%a.col(4))","(+3*a.col(11)+3*a.col(4)%a.col(10)+3*a.col(5)%a.col(1)+3*a.col(8)%a.col(8)+3*a.col(9)%a.col(0))"
    ),c("()","(+3*a.col(13))","(+3*a.col(0)%a.col(13)+3*a.col(8)%a.col(2)+6*a.col(1)%a.col(10)+6*a.col(10)%a.col(1)+3*a.col(2)%a.col(8)+3*a.col(13)%a.col(0))","(+3*a.col(11))","(+3*a.col(4)%a.col(11)+3*a.col(5)%a.col(10)+6*a.col(8)%a.col(9)+6*a.col(9)%a.col(8)+3*a.col(10)%a.col(5)+3*a.col(11)%a.col(4))","(+3*a.col(0)%a.col(11)+3*a.col(8)%a.col(10)+6*a.col(1)%a.col(9)+6*a.col(10)%a.col(8)+3*a.col(2)%a.col(5)+3*a.col(13)%a.col(4))","()","(+1*a.col(3))","(+1*a.col(0)%a.col(3)+3*a.col(1)%a.col(2)+3*a.col(2)%a.col(1)+1*a.col(3)%a.col(0))","(+1*a.col(13))","(+1*a.col(4)%a.col(13)+3*a.col(8)%a.col(10)+3*a.col(10)%a.col(8)+1*a.col(13)%a.col(4))","(+1*a.col(4)%a.col(3)+3*a.col(8)%a.col(2)+3*a.col(10)%a.col(1)+1*a.col(13)%a.col(0))","()","(+3*a.col(10))","(+3*a.col(13)+3*a.col(0)%a.col(10)+3*a.col(8)%a.col(1)+3*a.col(1)%a.col(8)+3*a.col(10)%a.col(0))","(+3*a.col(9))","(+3*a.col(12)+3*a.col(4)%a.col(9)+3*a.col(5)%a.col(8)+3*a.col(8)%a.col(5)+3*a.col(9)%a.col(4))","(+3*a.col(11)+3*a.col(4)%a.col(10)+3*a.col(5)%a.col(1)+3*a.col(8)%a.col(8)+3*a.col(9)%a.col(0))","()","(+1.5*a.col(2))","(+1.5*a.col(3)+1.5*a.col(0)%a.col(2)+3*a.col(1)%a.col(1)+1.5*a.col(2)%a.col(0))","(+1.5*a.col(10))","(+1.5*a.col(11)+1.5*a.col(4)%a.col(10)+3*a.col(8)%a.col(8)+1.5*a.col(10)%a.col(4))","(+1.5*a.col(13)+1.5*a.col(4)%a.col(2)+3*a.col(8)%a.col(1)+1.5*a.col(10)%a.col(0))","()","(+1.5*a.col(2))","(+1.5*a.col(3)+1.5*a.col(0)%a.col(2)+3*a.col(1)%a.col(1)+1.5*a.col(2)%a.col(0))","(+1.5*a.col(10))","(+1.5*a.col(11)+1.5*a.col(4)%a.col(10)+3*a.col(8)%a.col(8)+1.5*a.col(10)%a.col(4))","(+1.5*a.col(13)+1.5*a.col(4)%a.col(2)+3*a.col(8)%a.col(1)+1.5*a.col(10)%a.col(0))","()","()","()","()","()","()"
    ))

    MAT[1,1]='(1+0*a.col(0))'           # ONES(N2)
    MAT[5,7]='(1+0*a.col(0))'
    MAT[2,13]='(1+0*a.col(0))'
    MAT[9,19]='(1+0*a.col(0))'
    MAT[9,25]='(1+0*a.col(0))'
    MAT[6,31]='(1+0*a.col(0))'

   # HOMOGENEITY -----------------------------------------------------------------

   func.list.timehomo=func.list*0
   for(i in which(func.list==1))
   {
     # which expressions vary over time
      result=eval(body(namess[i]))
      func.list.timehomo[i]=2-(sum(diff(result)==0)==(length(result)-1))
   }
   if(any(func.list.timehomo==2)){homo=F}
   #func.list.timehomo[c(1+(0:5)*6)]=1    # Always set these to 1


   # BUILD ODE -------------------------------------------------------------------
   dims=rep('(',14)
   for(i in 1:14)
   {
     for(j in which(func.list==1))
     {
         if(MAT[i,j]!="()")
         {
           dims[i]=paste0(dims[i],'+(',body(namess[j])[2],')',c('*','%')[func.list.timehomo[j]],MAT[i,j])
         }
     }
     dims[i]=paste0(dims[i],')')
   }
     if(any(dims=='()')){dims[which(dims=='()')]='(1+0*a.col(0))'}
   for(i in 1:14)
   {
       dims[i]=paste0(paste0(paste0('        atemp.col(',i-1,')='),dims[i]),';')
   }

   # WRIGHT AND SOURCE -----------------------------------------------------------
     txt.full=paste(txtA,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],'\n',dims[6],'\n',dims[7],'\n',dims[8],'\n',dims[9],'\n',dims[10],'\n',dims[11],'\n',dims[12],'\n',dims[13],'\n',dims[14],txtB,txtC)
     type.sol ="                  GENERALIZED QUADRATIC DIFFUSON"
}

   #library(Rcpp)
   #library(RcppArmadillo)
   if(wrt)
   {
     write(txt.full,'BiGQD.mle.cpp')
   }
   stre="Compiling C++ code. Please wait."
   cat(stre, " \r")
   flush.console()
   sourceCpp(code=txt.full)
   cat('                                     ','\r')


#==============================================================================
#                           Interface Module
#==============================================================================
trim <- function (x) gsub("([[:space:]])", "", x)

namess4=c('a00','a10','a20','a01','a02','a11',
          'b00','b10','b20','b01','b02','b11',
          'c00','c10','c20','c01','c02','c11',
          'd00','d10','d20','d01','d02','d11',
          'e00','e10','e20','e01','e02','e11',
          'f00','f10','f20','f01','f02','f11')
indnames =rep(0,36)
for(i in 1:36)
{
  if(sum(obs==namess4[i]))
  {
    indnames[i]=TRUE
    namess4[i]=paste0(namess4[i],' : ',trim(body(namess4[i])[2]))
  }
}
namess4=matrix(namess4,length(namess4),1)
buffer0=c('================================================================')
buffer1=c('----------------------------------------------------------------')
buffer2=c('................................................................')
buffer3=c('...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ')
buffer4=c('_____________________ Drift Coefficients _______________________')
buffer5=c('___________________ Diffusion Coefficients _____________________')
buffer6=c('_____________________ Prior Distributions ______________________')
buffer7=c('__________________________ Model Info __________________________')

#Info1=data.frame(rbind(buffer0,matrix(names4[1:6,],6,1),buffer3,matrix(names4[1:6+6,],6,1),buffer3,matrix(names4[1:6+12,],6,1),buffer3,matrix(names4[1:6+18,],6,1),buffer3,matrix(names4[1:6+24,],6,1),buffer3,matrix(names4[1:6+30,],6,1),buffer2,prior.list,buffer),check.names=F)
#colnames(Info1)=type.sol
#buffer00=data.frame(buffer0,check.names=F)
#colnames(buffer00)=''
#print(buffer00 ,row.names = FALSE,right=F)
#print(Info1,row.names = FALSE,right=F)
# print(indnames)
Info=c(buffer0,type.sol,buffer0,buffer4,
       namess4[1:6][which(indnames[1:6]==T)],
       buffer3,
       namess4[7:12][which(indnames[7:12]==T)],
       buffer5,
       namess4[13:18][which(indnames[13:18]==T)],
       buffer3,
       namess4[19:24][which(indnames[19:24]==T)],
       buffer3,
       namess4[25:30][which(indnames[25:30]==T)],
       buffer3,
       namess4[31:36][which(indnames[31:36]==T)])
Info=data.frame(matrix(Info,length(Info),1))
colnames(Info)=''
if(print.output)
{
print(Info,row.names = FALSE,right=F)
}
############################################################################
############################################################################
############################################################################



    tme.eval = function(start_time)
    {
      start_time = as.POSIXct(start_time)
      dt = difftime(Sys.time(), start_time, units="secs")
      format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
     }
if(is.null(method))
{
  method = "Nelder-Mead"
}
if(is.null(control))
{
  control=list(fnscale=-1)
}
if(!is.null(control)) 
{
  if(any(names(control)=='fnscale'))
  {
    if(control$fnscale>=0)
    {
      stop(' ==============================================================================
           Incorrect input: control$fnscale must be negative!
      ==============================================================================')
    }      
    }
  if(!any(names(control)=='fnscale'))
  {
    control=c(control,fnscale=-1)
  }
  if(!any(names(control)=='maxit'))
  {
    control=c(control,maxit=5000)
  }
  }
     strts =matrix(0,nnn-1,2)
     strts[,1]=rtf[1]
     strts[,2]=rtf[2]
     like=function(pars)
     {
       rs=solver(X1[-nnn],X2[-nnn],X1[-1],X2[-1],c(0,pars),mesh,delt,nnn-1,T.seq[-nnn],strts,tro,secmom)
       sum(rs$like[-excl,1])
     }

     tme=Sys.time()
     result=optim(theta,like,control = control,method=method, hessian = T)
    tme=tme.eval(tme)



    actual.p=length(theta)

    model.inf=list(elapsed.time=tme,time.homogeneous=c('Yes','No')[2-homo],p=actual.p,N=length(X[,1])-length(excl)+1,Tag=Tag)
    Info2=c( buffer7,
             paste0("Time Homogeneous    : ",c('Yes','No')[2-homo]),
             paste0("Data Resolution     : ",c(paste0('Homogeneous: dt=',round(max(diff(T.seq)[-excl]),4)),paste0('Variable: min(dt)=',round(min(diff(T.seq)[-excl]),4),', max(dt)=',round(max(diff(T.seq)[-excl]),4)))[2-homo.res]),
             paste0("# Removed Transits. : ",c("None",length(excl))[2-is.null(exclude)]),
             paste0("Density approx.     : ",sol.state),
             paste0('Elapsed time        : ',tme),
             buffer3,
             paste0("dim(theta)          : ",round(actual.p,3)),
             buffer1)
    Info2=data.frame(matrix(Info2,length(Info2),1))
    colnames(Info2)=''
    if(print.output)
    {
    print(Info2,row.names = FALSE,right=F)
    }
    ret=list(opt=result,elapsed.time=tme,model.info=model.inf)
    class(ret)='GQD.mle'
    return(ret)
  }




