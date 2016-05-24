globalVariables('priors')
BiJGQD.mcmc=function(X,time,mesh=10,theta,sds,updates=10,burns=min(round(updates/2),25000),exclude=NULL,plot.chain=TRUE,RK.order=4,wrt=FALSE,Tag=NA,Dtype='Saddlepoint',Jdist='MVNormal',Jtype='Add',adapt=0,print.output=TRUE)
{
  Mstar3 =1
  recycle=FALSE
  rtf=runif(2)
  solver   =function(Xs, Xt, theta, N , delt , N2, tt  , P , alpha, lower , upper, tro  ){}
  rm(list =c('solver'))

  theta = theta+runif(length(theta),0.001,0.002)*sign(theta)
  T.seq=time
  Dtypes =c('Saddlepoint','Normal')
  Dindex = which(Dtypes==Dtype)

  JDtypes=c('Normal','MVNormal')
  JDindex = which(JDtypes==Jdist)

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
    ,'18. Density: Dtype has to be one of Saddlepoint or Normal.\n'
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
    ,'29. MCMC : Argument {burns} must be < {updates}.\n'
    ,'30. MCMC : Argument {updates} must be > 2.\n'
    ,'31. MCMC : length(theta)!=length(sds).\n'
    ,'32. Model: There has to be at least one model coefficient.\n'
    ,'33. Input: length(updates)!=1.\n'
    ,'34. Input: length(burns)!=1.\n'
    ,'35. Prior: priors(theta) must return a single value.\n'
    ,'36. Input: NAs not allowed.\n'
    ,'37. Input: length(Dtype)!=1.\n'
    ,'38. Input: NAs not allowed.\n'
    ,'39. Input: {Jdist} has to be of type "MVNormal".\n'
    ,'40. Input: {Jtype} has to be of type "Add" or "Mult".\n'
  )

   warntrue = rep(F,50)
    check.thetas = function(theta,tt)
  {
    t=tt
    theta = theta+runif(length(theta),0.001,0.002)*sign(theta)
    #namess=c('a00','a10','a20','a01','a02','a11',
    #         'b00','b10','b20','b01','b02','b11',
    #         'c00','c10','c20','c01','c02','c11',
    #         'd00','d10','d20','d01','d02','d11',
    #         'e00','e10','e20','e01','e02','e11',
    #         'f00','f10','f20','f01','f02','f11',
    #         'Lam00','Lamy0','Nmu11','Nsig11','Nmu21','Nsig21','Nmu12','Nsig12','Nmu22','Nsig22','Lam10','Lamx2','Lamy1','Lamy2','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')
    namess=c('a00','a10','a20','a01','a02','a11',
             'b00','b10','b20','b01','b02','b11',
             'c00','c10','c20','c01','c02','c11',
             'd00','d10','d20','d01','d02','d11',
             'e00','e10','e20','e01','e02','e11',
             'f00','f10','f20','f01','f02','f11',
             'Lam00','Lam10','Lam01','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')
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
        theta[j] = theta[j]+runif(1,0.1,1)
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
             'f00','f10','f20','f01','f02','f11',
             'Lam00','Lam10','Lam01','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')
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
   if(missing(sds))                                              {warntrue[4]=T}
   if(!is.matrix(X))                                             {warntrue[5]=T}
   if(!is.vector(time))                                          {warntrue[6]=T}
   # Check model parameters:
   if(check.thetas2(theta)!=0)                                   {warntrue[7]=T}
   if(!warntrue[7]){if(any(check.thetas(theta,T.seq)==0))        {warntrue[8]=T}}

   # Check input length:
   if(dim(X)[1]<10)                                              {warntrue[9]=T}
   if(length(time)<10)                                          {warntrue[10]=T}
   if(length(mesh)!=1)                                        {warntrue[14]=T}
   if(length(RK.order)!=1)                                      {warntrue[17]=T}
   if(length(updates)!=1)                                     {warntrue[33]=T}
   if(length(burns)!=1)                                         {warntrue[34]=T}
   if(length(Dtype)!=1)                                         {warntrue[37]=T}


   # Check density approx parameters:
   if(sum(Dindex)==0)                                          {warntrue[18] =T}


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
   if(!any(warntrue[c(33,34)])){if(burns>updates)             {warntrue[29]=T}}
   if(!warntrue[33]){if(updates<2)                            {warntrue[30]=T}}
   if(length(theta)!=length(sds))                             {warntrue[31]=T}
   if(any(is.na(X))||any(is.na(time)))                        {warntrue[36]=T}
   if(sum(JDindex)==0)                                        {warntrue[39] =TRUE}
   if((Jtype!='Add')&&(Jtype!='Mult'))                        {warntrue[40] =TRUE}
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
   #                           Prior distribution Module
   #==============================================================================
    if(sum(objects(pos=1)=='priors')==1)
    {
       pp=function(theta){}
       body(pp)=parse(text =body(priors)[2])
       prior.list=paste0('d(theta)',':',paste0(body(priors)[2]))
       if(length(priors(theta))!=1){stop(" ==============================================================================
   Incorrect input: Prior distribution must return a single value only!
   ==============================================================================");}
    }
    if(sum(objects(pos=1)=='priors')==0)
    {
      prior.list=paste0('d(theta)',':',' None.')
      pp=function(theta){1}
    }
    if(length(pp(theta))!=1){stop("Prior density must return only a single value!")}
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
    state1=F
    state3.0=F
    state3.1=F
    state3.2=F
    state4=T
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
    tro = 14


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
        mat atemp(N2,34);'
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
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N,double delt,int N2,vec tt,mat starts,int tro,int secmom,vec seq1,vec seq2)
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
    for (int i = 1; i <= 14; i++)
    {
     x0.col(i-1)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
     x0.col(i-1+14)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
    }

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
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N, mat delt,int N2,vec d,mat starts,int tro,int secmom,vec seq1,vec seq2)
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
    for (int i = 1; i <= 14; i++)
    {
     x0.col(i-1)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
     x0.col(i-1+14)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
    }
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
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N,double delt,int N2,vec tt,mat starts,int tro,int secmom,vec seq1,vec seq2)
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
  double whch =0;

  for (int i = 1; i <= 14; i++)
  {
   x0.col(i-1)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
   x0.col(i-1+14)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
  }

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
List  solver(vec Xs,vec Ys,vec Xt,vec Yt,vec theta,int N,double delt,int N2,vec tt,mat starts,int tro,int secmom,vec seq1,vec seq2)
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
  double whch =0;
  x0.fill(0);
  for (int i = 1; i <= 14; i++)
  {
     x0.col(i-1)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
     x0.col(i-1+14)=pow(Xs,seq1[i-1])%pow(Ys,seq2[i-1]);
  }

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


if(state4)
{
 if(Dtype=='Saddlepoint')
 {
 txtC='
 vec probs=exp(-x0.col(28)-x0.col(29)-x0.col(30)-x0.col(31)-x0.col(32)-x0.col(33));

   vec m00 = (1+0*x0.col(0));
   vec m10 = (x0.col(0) -probs%x0.col(14))/(1-probs);
   vec m20 = (x0.col(1) -probs%x0.col(15))/(1-probs);
   vec m30 = (x0.col(2) -probs%x0.col(16))/(1-probs);
   vec m40 = (x0.col(3) -probs%x0.col(17))/(1-probs);
   vec m01 = (x0.col(4) -probs%x0.col(18))/(1-probs);
   vec m02 = (x0.col(5) -probs%x0.col(19))/(1-probs);
   vec m03 = (x0.col(6) -probs%x0.col(20))/(1-probs);
   vec m04 = (x0.col(7) -probs%x0.col(21))/(1-probs);
   vec m11 = (x0.col(8) -probs%x0.col(22))/(1-probs);
   vec m12 = (x0.col(9) -probs%x0.col(23))/(1-probs);
   vec m21 = (x0.col(10)-probs%x0.col(24))/(1-probs);
   vec m22 = (x0.col(11)-probs%x0.col(25))/(1-probs);
   vec m13 = (x0.col(12)-probs%x0.col(26))/(1-probs);
   vec m31 = (x0.col(13)-probs%x0.col(27))/(1-probs);

   vec mm00 = (1+0*x0.col(0));
   vec mm10 = x0.col(14);
   vec mm20 = x0.col(15);
   vec mm30 = x0.col(16);
   vec mm40 = x0.col(17);
   vec mm01 = x0.col(18);
   vec mm02 = x0.col(19);
   vec mm03 = x0.col(20);
   vec mm04 = x0.col(21);
   vec mm11 = x0.col(22);
   vec mm12 = x0.col(23);
   vec mm21 = x0.col(24);
   vec mm22 = x0.col(25);
   vec mm13 = x0.col(26);
   vec mm31 = x0.col(27);


  vec k10 = m10;
  vec k20 = m20-pow(m10,2);
  vec k30 = m30-3*m20%m10+2*pow(m10,3);
  vec k40 = m40 -4*m30%m10-3*pow(m20,2)+12*m20%pow(m10,2)-6*pow(m10,4);
  vec k01 = m01;
  vec k02 = m02-  pow(m01,2) ;
  vec k03 = m03-3*m02%m01+2*pow(m01,3) ;
  vec k04 = m04-4*m03%m01-3*pow(m02,2)+12*m02%pow(m01,2)-6*pow(m01,4);
  vec k11 = m11-m10%m01;
  vec k21 = m21-2*m11%m10-m20%m01+2*pow(m10,2)%m01;
  vec k12 = m12-2*m11%m01-m02%m10+2*pow(m01,2)%m10;
  vec k22 = m22-2*m21%m01-2*m12%m10-m20%m02-2*pow(m11,2)+8*m11%m01%m10+2*m02%pow(m10,2)+2*m20%pow(m01,2)-6*pow(m10,2)%pow(m01,2) ;
  vec k31 = m31-3*m21%m10-m30%m01-3*m20%m11+6*m11%pow(m10,2)+6*m20%m10%m01-6*pow(m10,3)%m01 ;
  vec k13 = m13-3*m12%m01-m03%m10-3*m02%m11+6*m11%pow(m01,2)+6*m02%m01%m10-6*pow(m01,3)%m10;

  vec kk10 = mm10;
  vec kk20 = mm20-pow(mm10,2);
  vec kk30 = mm30-3*mm20%mm10+2*pow(mm10,3);
  vec kk40 = mm40 -4*mm30%mm10-3*pow(mm20,2)+12*mm20%pow(mm10,2)-6*pow(mm10,4);
  vec kk01 = mm01;
  vec kk02 = mm02-  pow(mm01,2) ;
  vec kk03 = mm03-3*mm02%mm01+2*pow(mm01,3) ;
  vec kk04 = mm04-4*mm03%mm01-3*pow(mm02,2)+12*mm02%pow(mm01,2)-6*pow(mm01,4);
  vec kk11 = mm11-mm10%mm01;
  vec kk21 = mm21-2*mm11%mm10-mm20%mm01+2*pow(mm10,2)%mm01;
  vec kk12 = mm12-2*mm11%mm01-mm02%mm10+2*pow(mm01,2)%mm10;
  vec kk22 = mm22-2*mm21%mm01-2*mm12%mm10-mm20%mm02-2*pow(mm11,2)+8*mm11%mm01%mm10+2*mm02%pow(mm10,2)+2*mm20%pow(mm01,2)-6*pow(mm10,2)%pow(mm01,2) ;
  vec kk31 = mm31-3*mm21%mm10-mm30%mm01-3*mm20%mm11+6*mm11%pow(mm10,2)+6*mm20%mm10%mm01-6*pow(mm10,3)%mm01 ;
  vec kk13 = mm13-3*mm12%mm01-mm03%mm10-3*mm02%mm11+6*mm11%pow(mm01,2)+6*mm02%mm01%mm10-6*pow(mm01,3)%mm10;


  vec a(N2);
  vec b(N2);
  vec abser(N2);
  abser=0.1+abser;
  a.ones();
  b.ones();
  vec det=(k10%k01-k11%k11);
  a=-(Xt-k10)%k20/det+(Yt-k01)%k11/det;
  b=+(Xt-k10)%k11/det-(Yt-k01)%k02/det;
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
  while((max(abser)>0.001)&&(ind<1500))
  {
  gg=k10+k20%a+(1.0/2.0)*k30%a%a+(1.0/6.0)*k40%a%a%a +k11%b +(1.0/2.0)*k12%b%b+k21%a%b+(1.0/6.0)*b%b%b%k13+(1.0/2.0)*a%a%b%k31+(1.0/2.0)*a%b%b%k22-Xt;
  hh=k01+k02%b+(1.0/2.0)*k03%b%b+(1.0/6.0)*k04%b%b%b +k11%a +k12%a%b+(1.0/2.0)*k21%a%a+(1.0/2.0)*a%b%b%k13+(1.0/6.0)*a%a%a%k31+(1.0/2.0)*a%a%b%k22-Yt;

  gg1=k20+k30%a+(1.0/2.0)*k40%a%a+k21%b+a%b%k31+(1.0/2.0)*b%b%k22;
  gg2=k11 +k12%b+k21%a+(1.0/2.0)*b%b%k13+(1.0/2.0)*a%a%k31+a%b%k22;

  hh1=k11 +k12%b+k21%a+(1.0/2.0)*b%b%k13+(1.0/2.0)*a%a%k31+a%b%k22;
  hh2=k02+k03%b+(1.0/2.0)*k04%b%b +k12%a+a%b%k13+(1.0/2.0)*a%a%k22;

  anew  =a+(hh%gg2-gg%hh2)/(gg1%hh2-gg2%hh1);
  bnew  =b-(hh%gg1-gg%hh1)/(gg1%hh2-gg2%hh1);
  abser =(pow(anew-a,2)+pow(bnew-b,2));
  a=anew;
  b=bnew;
  ind=ind+1;
  }

  vec dens2=-log(2*3.141592653589793)-0.5*log(gg1%hh2-gg2%hh1)+(k10%a+k01%b+(1.0/2.0)*k20%a%a+(1.0/2.0)*k02%b%b+(1.0/6.0)*k30%a%a%a+(1.0/6.0)*k03%b%b%b+(1.0/24.0)*k40%a%a%a%a+(1.0/24.0)*k04%b%b%b%b+k11%a%b+(1.0/2.0)*k12%a%b%b+(1.0/2.0)*k21%a%a%b+(1.0/6.0)*a%b%b%b%k13+(1.0/6.0)*a%a%a%b%k31+(1.0/4.0)*a%a%b%b%k22-a%Xt-b%Yt);

  abser=0.1+0*abser;
  a.ones();
  b.ones();
  det=(kk10%kk01-kk11%kk11);
  a=-(Xt-kk10)%kk20/det+(Yt-kk01)%kk11/det;
  b=+(Xt-kk10)%kk11/det-(Yt-kk01)%kk02/det;
  ind=0;
  while((max(abser)>0.001)&&(ind<1500))
  {
  gg=kk10+kk20%a+(1.0/2.0)*kk30%a%a+(1.0/6.0)*kk40%a%a%a +kk11%b +(1.0/2.0)*kk12%b%b+kk21%a%b+(1.0/6.0)*b%b%b%kk13+(1.0/2.0)*a%a%b%kk31+(1.0/2.0)*a%b%b%kk22-Xt;
  hh=kk01+kk02%b+(1.0/2.0)*kk03%b%b+(1.0/6.0)*kk04%b%b%b +kk11%a +kk12%a%b+(1.0/2.0)*kk21%a%a+(1.0/2.0)*a%b%b%kk13+(1.0/6.0)*a%a%a%kk31+(1.0/2.0)*a%a%b%kk22-Yt;

  gg1=kk20+kk30%a+(1.0/2.0)*kk40%a%a+kk21%b+a%b%kk31+(1.0/2.0)*b%b%kk22;
  gg2=kk11 +kk12%b+kk21%a+(1.0/2.0)*b%b%kk13+(1.0/2.0)*a%a%kk31+a%b%kk22;

  hh1=kk11 +kk12%b+kk21%a+(1.0/2.0)*b%b%kk13+(1.0/2.0)*a%a%kk31+a%b%kk22;
  hh2=kk02+kk03%b+(1.0/2.0)*kk04%b%b +kk12%a+a%b%kk13+(1.0/2.0)*a%a%kk22;

  anew  =a+(hh%gg2-gg%hh2)/(gg1%hh2-gg2%hh1);
  bnew  =b-(hh%gg1-gg%hh1)/(gg1%hh2-gg2%hh1);
  abser =(pow(anew-a,2)+pow(bnew-b,2));
  a=anew;
  b=bnew;
  ind=ind+1;
  }
  vec dens1=-log(2*3.141592653589793)-0.5*log(gg1%hh2-gg2%hh1)+(kk10%a+kk01%b+(1.0/2.0)*kk20%a%a+(1.0/2.0)*kk02%b%b+(1.0/6.0)*kk30%a%a%a+(1.0/6.0)*kk03%b%b%b+(1.0/24.0)*kk40%a%a%a%a+(1.0/24.0)*kk04%b%b%b%b+kk11%a%b+(1.0/2.0)*kk12%a%b%b+(1.0/2.0)*kk21%a%a%b+(1.0/6.0)*a%b%b%b%kk13+(1.0/6.0)*a%a%a%b%kk31+(1.0/4.0)*a%a%b%b%kk22-a%Xt-b%Yt);

  resss.col(1)=a;
  resss.col(2)=b;
  resss.col(0)=log(exp(dens1)%probs+exp(dens2)%(1-probs));
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  ret["probs"] = probs;
  return(ret);  }'
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
  vec probs=exp(-x0.col(28)-x0.col(29)-x0.col(30)-x0.col(31)-x0.col(32)-x0.col(33));

   vec m00 = (1+0*x0.col(0));
   vec m10 = (x0.col(0) -probs%x0.col(14))/(1-probs);
   vec m20 = (x0.col(1) -probs%x0.col(15))/(1-probs);
   vec m30 = (x0.col(2) -probs%x0.col(16))/(1-probs);
   vec m40 = (x0.col(3) -probs%x0.col(17))/(1-probs);
   vec m01 = (x0.col(4) -probs%x0.col(18))/(1-probs);
   vec m02 = (x0.col(5) -probs%x0.col(19))/(1-probs);
   vec m03 = (x0.col(6) -probs%x0.col(20))/(1-probs);
   vec m04 = (x0.col(7) -probs%x0.col(21))/(1-probs);
   vec m11 = (x0.col(8) -probs%x0.col(22))/(1-probs);
   vec m12 = (x0.col(9) -probs%x0.col(23))/(1-probs);
   vec m21 = (x0.col(10)-probs%x0.col(24))/(1-probs);
   vec m22 = (x0.col(11)-probs%x0.col(25))/(1-probs);
   vec m13 = (x0.col(12)-probs%x0.col(26))/(1-probs);
   vec m31 = (x0.col(13)-probs%x0.col(27))/(1-probs);

   vec mm00 = (1+0*x0.col(0));
   vec mm10 = x0.col(14);
   vec mm20 = x0.col(15);
   vec mm30 = x0.col(16);
   vec mm40 = x0.col(17);
   vec mm01 = x0.col(18);
   vec mm02 = x0.col(19);
   vec mm03 = x0.col(20);
   vec mm04 = x0.col(21);
   vec mm11 = x0.col(22);
   vec mm12 = x0.col(23);
   vec mm21 = x0.col(24);
   vec mm22 = x0.col(25);
   vec mm13 = x0.col(26);
   vec mm31 = x0.col(27);

  vec k10 = m10;
  vec k20 = m20-pow(m10,2);
  vec k30 = m30-3*m20%m10+2*pow(m10,3);
  vec k40 = m40 -4*m30%m10-3*pow(m20,2)+12*m20%pow(m10,2)-6*pow(m10,4);
  vec k01 = m01;
  vec k02 = m02-  pow(m01,2) ;
  vec k03 = m03-3*m02%m01+2*pow(m01,3) ;
  vec k04 = m04-4*m03%m01-3*pow(m02,2)+12*m02%pow(m01,2)-6*pow(m01,4);
  vec k11 = m11-m10%m01;
  vec k21 = m21-2*m11%m10-m20%m01+2*pow(m10,2)%m01;
  vec k12 = m12-2*m11%m01-m02%m10+2*pow(m01,2)%m10;
  vec k22 = m22-2*m21%m01-2*m12%m10-m20%m02-2*pow(m11,2)+8*m11%m01%m10+2*m02%pow(m10,2)+2*m20%pow(m01,2)-6*pow(m10,2)%pow(m01,2) ;
  vec k31 = m31-3*m21%m10-m30%m01-3*m20%m11+6*m11%pow(m10,2)+6*m20%m10%m01-6*pow(m10,3)%m01 ;
  vec k13 = m13-3*m12%m01-m03%m10-3*m02%m11+6*m11%pow(m01,2)+6*m02%m01%m10-6*pow(m01,3)%m10;

  vec kk10 = mm10;
  vec kk20 = mm20-pow(mm10,2);
  vec kk30 = mm30-3*mm20%mm10+2*pow(mm10,3);
  vec kk40 = mm40 -4*mm30%mm10-3*pow(mm20,2)+12*mm20%pow(mm10,2)-6*pow(mm10,4);
  vec kk01 = mm01;
  vec kk02 = mm02-  pow(mm01,2) ;
  vec kk03 = mm03-3*mm02%mm01+2*pow(mm01,3) ;
  vec kk04 = mm04-4*mm03%mm01-3*pow(mm02,2)+12*mm02%pow(mm01,2)-6*pow(mm01,4);
  vec kk11 = mm11-mm10%mm01;
  vec kk21 = mm21-2*mm11%mm10-mm20%mm01+2*pow(mm10,2)%mm01;
  vec kk12 = mm12-2*mm11%mm01-mm02%mm10+2*pow(mm01,2)%mm10;
  vec kk22 = mm22-2*mm21%mm01-2*mm12%mm10-mm20%mm02-2*pow(mm11,2)+8*mm11%mm01%mm10+2*mm02%pow(mm10,2)+2*mm20%pow(mm01,2)-6*pow(mm10,2)%pow(mm01,2) ;
  vec kk31 = mm31-3*mm21%mm10-mm30%mm01-3*mm20%mm11+6*mm11%pow(mm10,2)+6*mm20%mm10%mm01-6*pow(mm10,3)%mm01 ;
  vec kk13 = mm13-3*mm12%mm01-mm03%mm10-3*mm02%mm11+6*mm11%pow(mm01,2)+6*mm02%mm01%mm10-6*pow(mm01,3)%mm10;


  vec det=(kk20%kk02-kk11%kk11);
  vec dens1=-log(2*3.141592653589793)-0.5*log(abs(det))-0.5*((Xt-kk10)%(Xt-kk10)%kk02/det-(Xt-kk10)%(Yt-kk01)%kk11/det-(Xt-kk10)%(Yt-kk01)%kk11/det+(Yt-kk01)%(Yt-kk01)%kk20/det);
  det=(k20%k02-k11%k11);
  vec dens2=-log(2*3.141592653589793)-0.5*log(abs(det))-0.5*((Xt-k10)%(Xt-k10)%k02/det-(Xt-k10)%(Yt-k01)%k11/det-(Xt-k10)%(Yt-k01)%k11/det+(Yt-k01)%(Yt-k01)%k20/det);

  resss.col(0)=log(exp(dens1)%probs+exp(dens2)%(1-probs));
  List ret;
  ret["like"] = resss;
  ret["max"] = whch;
  ret["probs"] = probs;
  return(ret);
}'
}
}

   #==============================================================================
   #                   Generate TYPE of Solution
   #==============================================================================

if(state4)
{
   # DATA RESOLUTION -------------------------------------------------------------
   if(!homo.res)
   {
     delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,
                diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,
                diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,
                diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,
                diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh)

   }

  t=seq(0,100,1/100)
  tro=34

  if(Jdist=='Normal')
  {
    Mstar1 =
   c("+1*mm1*m00"
   ,"+2*mm1*m10+1*mm2"
   ,"+3*mm1*m20+3*mm2*m10+1*mm3"
   ,"+4*mm1*m30+6*mm2*m20+4*mm3*m10+1*mm4"
   ,"+1*oo1*m00"
   ,"+2*oo1*m01+1*oo2"
   ,"+3*oo1*m02+3*oo2*m01+1*oo3"
   ,"+4*oo1*m03+6*oo2*m02+4*oo3*m01+1*oo4"
   ,"mm1*m01+oo1*m10+mm1*oo1"
   ,"mm1*m02+2*oo1*m11+2*mm1*oo1*m01+oo2*m10+mm1*oo2"
   ,"2*mm1*m11+mm2*m01+oo1*m20+2*mm1*oo1*m10+mm2*oo1"
   ,"2*mm1*m12+mm2*m02+2*oo1*m21+4*mm1*oo1*m11+2*mm2*oo1*m01+oo2*m20+2*mm1*oo2*m10+mm2*oo2"
   ,"mm1*m03+3*oo1*m12+3*mm1*oo1*m02+3*oo2*m11+3*mm1*oo2*m01+oo3*m10+mm1*oo3"
   ,"3*mm1*m21+3*mm2*m11+mm3*m01+oo1*m30+3*mm1*oo1*m20+3*mm2*oo1*m10+mm3*oo1")

    Mstar2 =
   c("+1*nn1*m00"
   ,"+2*nn1*m10+1*nn2"
   ,"+3*nn1*m20+3*nn2*m10+1*nn3"
   ,"+4*nn1*m30+6*nn2*m20+4*nn3*m10+1*nn4"
   ,"+1*pp1*m00"
   ,"+2*pp1*m01+1*pp2"
   ,"+3*pp1*m02+3*pp2*m01+1*pp3"
   ,"+4*pp1*m03+6*pp2*m02+4*pp3*m01+1*pp4"
   ,"nn1*m01+pp1*m10+nn1*pp1"
   ,"nn1*m02+2*pp1*m11+2*nn1*pp1*m01+pp2*m10+nn1*pp2"
   ,"2*nn1*m11+nn2*m01+pp1*m20+2*nn1*pp1*m10+nn2*pp1"
   ,"2*nn1*m12+nn2*m02+2*pp1*m21+4*nn1*pp1*m11+2*nn2*pp1*m01+pp2*m20+2*nn1*pp2*m10+nn2*pp2"
   ,"nn1*m03+3*pp1*m12+3*nn1*pp1*m02+3*pp2*m11+3*nn1*pp2*m01+pp3*m10+nn1*pp3"
   ,"3*nn1*m21+3*nn2*m11+nn3*m01+pp1*m30+3*nn1*pp1*m20+3*nn2*pp1*m10+nn3*pp1")

      Mstar1x =
   c("+1*mm1*m10"
   ,"+2*mm1*m20+1*mm2*m10"
   ,"+3*mm1*m30+3*mm2*m20+1*mm3*m10"
   ,"+4*mm1*m40+6*mm2*m30+4*mm3*m20+1*mm4*m10"
   ,"+1*oo1*m10"
   ,"+2*oo1*m11+1*oo2*m10"
   ,"+3*oo1*m12+3*oo2*m11+1*oo3*m10"
   ,"+4*oo1*m13+6*oo2*m12+4*oo3*m11+1*oo4*m10"
   ,"mm1*m11+oo1*m20+mm1*oo1*m10"
   ,"mm1*m12+2*oo1*m21+2*mm1*oo1*m11+oo2*m20+mm1*oo2*m10"
   ,"2*mm1*m21+mm2*m11+oo1*m30+2*mm1*oo1*m20+mm2*oo1*m10"
   ,"2*mm1*m22+mm2*m12+2*oo1*m31+4*mm1*oo1*m21+2*mm2*oo1*m11+oo2*m30+2*mm1*oo2*m20+mm2*oo2*m10"
   ,"mm1*m13+3*oo1*m22+3*mm1*oo1*m12+3*oo2*m21+3*mm1*oo2*m11+oo3*m20+mm1*oo3*m10"
   ,"3*mm1*m31+3*mm2*m21+mm3*m11+oo1*m40+3*mm1*oo1*m30+3*mm2*oo1*m20+mm3*oo1*m10")

    Mstar1y =
   c("+1*mm1*m01"
   ,"+2*mm1*m11+1*mm2*m01"
   ,"+3*mm1*m21+3*mm2*m11+1*mm3*m01"
   ,"+4*mm1*m31+6*mm2*m21+4*mm3*m11+1*mm4*m01"
   ,"+1*oo1*m01"
   ,"+2*oo1*m02+1*oo2*m01"
   ,"+3*oo1*m03+3*oo2*m02+1*oo3*m01"
   ,"+4*oo1*m04+6*oo2*m03+4*oo3*m02+1*oo4*m01"
   ,"mm1*m02+oo1*m11+mm1*oo1*m01"
   ,"mm1*m03+2*oo1*m11+2*mm1*oo1*m02+oo2*m11+mm1*oo2*m01"
   ,"2*mm1*m12+mm2*m02+oo1*m21+2*mm1*oo1*m11+mm2*oo1*m01"
   ,"2*mm1*m13+mm2*m03+2*oo1*m22+4*mm1*oo1*m12+2*mm2*oo1*m02+oo2*m21+2*mm1*oo2*m11+mm2*oo2*m01"
   ,"mm1*m04+3*oo1*m13+3*mm1*oo1*m03+3*oo2*m11+3*mm1*oo2*m02+oo3*m11+mm1*oo3*m01"
   ,"3*mm1*m22+3*mm2*m12+mm3*m02+oo1*m31+3*mm1*oo1*m21+3*mm2*oo1*m11+mm3*oo1*m02")

    Mstar2x =
   c("+1*nn1*m10"
   ,"+2*nn1*m20+1*nn2*m10"
   ,"+3*nn1*m30+3*nn2*m20+1*nn3*m10"
   ,"+4*nn1*m40+6*nn2*m30+4*nn3*m20+1*nn4*m10"
   ,"+1*pp1*m10"
   ,"+2*pp1*m11+1*pp2*m10"
   ,"+3*pp1*m12+3*pp2*m11+1*pp3*m10"
   ,"+4*pp1*m13+6*pp2*m12+4*pp3*m11+1*pp4*m10"
   ,"nn1*m11+pp1*m20+nn1*pp1*m10"
   ,"nn1*m12+2*pp1*m21+2*nn1*pp1*m11+pp2*m20+nn1*pp2*m10"
   ,"2*nn1*m21+nn2*m11+pp1*m30+2*nn1*pp1*m20+nn2*pp1*m10"
   ,"2*nn1*m22+nn2*m12+2*pp1*m31+4*nn1*pp1*m21+2*nn2*pp1*m11+pp2*m30+2*nn1*pp2*m20+nn2*pp2*m10"
   ,"nn1*m13+3*pp1*m22+3*nn1*pp1*m12+3*pp2*m21+3*nn1*pp2*m11+pp3*m20+nn1*pp3*m10"
   ,"3*nn1*m31+3*nn2*m21+nn3*m11+pp1*m40+3*nn1*pp1*m30+3*nn2*pp1*m20+nn3*pp1*m10")


    Mstar2y =
   c("+1*nn1*m01"
   ,"+2*nn1*m11+1*nn2*m01"
   ,"+3*nn1*m21+3*nn2*m11+1*nn3*m01"
   ,"+4*nn1*m31+6*nn2*m21+4*nn3*m11+1*nn4*m01"
   ,"+1*pp1*m01"
   ,"+2*pp1*m02+1*pp2*m01"
   ,"+3*pp1*m03+3*pp2*m02+1*pp3*m01"
   ,"+4*pp1*m04+6*pp2*m03+4*pp3*m02+1*pp4*m01"
   ,"nn1*m02+pp1*m11+nn1*pp1*m01"
   ,"nn1*m03+2*pp1*m12+2*nn1*pp1*m02+pp2*m11+nn1*pp2*m01"
   ,"2*nn1*m12+nn2*m02+pp1*m21+2*nn1*pp1*m11+nn2*pp1*m01"
   ,"2*nn1*m13+nn2*m03+2*pp1*m22+4*nn1*pp1*m12+2*nn2*pp1*m02+pp2*m21+2*nn1*pp2*m11+nn2*pp2*m01"
   ,"nn1*m04+3*pp1*m13+3*nn1*pp1*m03+3*pp2*m12+3*nn1*pp2*m02+pp3*m11+nn1*pp3*m01"
   ,"3*nn1*m22+3*nn2*m12+nn3*m02+pp1*m31+3*nn1*pp1*m21+3*nn2*pp1*m11+nn3*pp1*m01")
  }

  if(Jdist=='MVNormal')
  {
    if(Jtype=='Add')
   {
     Mstar1 =
  c("+1*mv10*m00"
   ,"+2*mv10*m10+1*mv20"
   ,"+3*mv10*m20+3*mv20*m10+1*mv30"
   ,"+4*mv10*m30+6*mv20*m20+4*mv30*m10+1*mv40"
   ,"+1*mv01*m00"
   ,"+2*mv01*m01+1*mv02"
   ,"+3*mv01*m02+3*mv02*m01+1*mv03"
   ,"+4*mv01*m03+6*mv02*m02+4*mv03*m01+1*mv04"
   ,"mv10*m01+mv01*m10+mv11"
   ,"mv10*m02+2*mv01*m11+2*mv11*m01+mv02*m10+mv12"
   ,"2*mv10*m11+mv20*m01+mv01*m20+2*mv11*m10+mv21"
   ,"2*mv10*m12+mv20*m02+2*mv01*m21+4*mv11*m11+2*mv21*m01+mv02*m20+2*mv12*m10+mv22"
   ,"mv10*m03+3*mv01*m12+3*mv11*m02+3*mv02*m11+3*mv12*m01+mv03*m10+mv13"
   ,"3*mv10*m21+3*mv20*m11+mv30*m01+mv01*m30+3*mv11*m20+3*mv21*m10+mv31")
   }
   if(Jtype=='Mult')
   {
     Mstar1 =
  c("m10*(mv10)"
   ,"m20*(mv20+2*mv10)"
   ,"m30*(mv30+3*mv20+3*mv10)"
   ,"m40*(mv40+4*mv30+6*mv20+4*mv10)"
   ,"m01*(mv01)"
   ,"m02*(mv02+2*mv01)"
   ,"m03*(mv03+3*mv02+3*mv01)"
   ,"m04*(mv04+4*mv03+6*mv02+4*mv01)"
   ,"m11*(mv11+mv01+mv10)"
   ,"m12*(mv12+mv02+2*mv11+2*mv01+mv10)"
   ,"m21*(mv21+2*mv11+mv01+mv20+2*mv10)"
   ,"m22*(mv22+2*mv12+mv02+2*mv21+4*mv11+2*mv01+mv20+2*mv10)"
   ,"m13*(mv13+mv03+3*mv12+3*mv02+3*mv11+3*mv01+mv10)"
   ,"m31*(mv31+3*mv21+3*mv11+mv01+mv30+3*mv20+3*mv10)")
   }

  }
MAT=rbind(
 c("1*m00","1*m10","1*m20","1*m01","1*m02","1*m11","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","")
,c("2*m10","2*m20","2*m30","2*m11","2*m12","2*m21","","","","","","","1*m00","1*m10","1*m20","1*m01","1*m02","1*m11","","","","","","","","","","","","","","","","","","")
,c("3*m20","3*m30","3*m40","3*m21","3*m22","3*m31","","","","","","","3*m10","3*m20","3*m30","3*m11","3*m12","3*m21","","","","","","","","","","","","","","","","","","")
,c("4*m30","4*m40","4*m50","4*m31","4*m32","4*m41","","","","","","","6*m20","6*m30","6*m40","6*m21","6*m22","6*m31","","","","","","","","","","","","","","","","","","")
,c("","","","","","","1*m00","1*m10","1*m20","1*m01","1*m02","1*m11","","","","","","","","","","","","","","","","","","","","","","","","")
,c("","","","","","","2*m01","2*m11","2*m21","2*m02","2*m03","2*m12","","","","","","","","","","","","","","","","","","","1*m00","1*m10","1*m20","1*m01","1*m02","1*m11")
,c("","","","","","","3*m02","3*m12","3*m22","3*m03","3*m04","3*m13","","","","","","","","","","","","","","","","","","","3*m01","3*m11","3*m21","3*m02","3*m03","3*m12")
,c("","","","","","","4*m03","4*m13","4*m23","4*m04","4*m05","4*m14","","","","","","","","","","","","","","","","","","","6*m02","6*m12","6*m22","6*m03","6*m04","6*m13")
,c("1*m01","1*m11","1*m21","1*m02","1*m03","1*m12","1*m10","1*m20","1*m30","1*m11","1*m12","1*m21","","","","","","","1*m00","1*m10","1*m20","1*m01","1*m02","1*m11","1*m00","1*m10","1*m20","1*m01","1*m02","1*m11","","","","","","")
,c("1*m02","1*m12","1*m22","1*m03","1*m04","1*m13","2*m11","2*m21","2*m31","2*m12","2*m13","2*m22","","","","","","","2*m01","2*m11","2*m21","2*m02","2*m03","2*m12","2*m01","2*m11","2*m21","2*m02","2*m03","2*m12","1*m10","1*m20","1*m30","1*m11","1*m12","1*m21")
,c("2*m11","2*m21","2*m31","2*m12","2*m13","2*m22","1*m20","1*m30","1*m40","1*m21","1*m22","1*m31","1*m01","1*m11","1*m21","1*m02","1*m03","1*m12","2*m10","2*m20","2*m30","2*m11","2*m12","2*m21","2*m10","2*m20","2*m30","2*m11","2*m12","2*m21","","","","","","")
,c("2*m12","2*m22","2*m32","2*m13","2*m14","2*m23","2*m21","2*m31","2*m41","2*m22","2*m23","2*m32","1*m02","1*m12","1*m22","1*m03","1*m04","1*m13","4*m11","4*m21","4*m31","4*m12","4*m13","4*m22","4*m11","4*m21","4*m31","4*m12","4*m13","4*m22","1*m20","1*m30","1*m40","1*m21","1*m22","1*m31")
,c("1*m03","1*m13","1*m23","1*m04","1*m05","1*m14","3*m12","3*m22","3*m32","3*m13","3*m14","3*m23","","","","","","","3*m02","3*m12","3*m22","3*m03","3*m04","3*m13","3*m02","3*m12","3*m22","3*m03","3*m04","3*m13","3*m11","3*m21","3*m31","3*m12","3*m13","3*m22")
,c("3*m21","3*m31","3*m41","3*m22","3*m23","3*m32","1*m30","1*m40","1*m50","1*m31","1*m32","1*m41","3*m11","3*m21","3*m31","3*m12","3*m13","3*m22","3*m20","3*m30","3*m40","3*m21","3*m22","3*m31","3*m20","3*m30","3*m40","3*m21","3*m22","3*m31","","","","","","")
)

MAT2=rbind(
 c("1*mm00","1*mm10","1*mm20","1*mm01","1*mm02","1*mm11","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","")
,c("2*mm10","2*mm20","2*mm30","2*mm11","2*mm12","2*mm21","","","","","","","1*mm00","1*mm10","1*mm20","1*mm01","1*mm02","1*mm11","","","","","","","","","","","","","","","","","","")
,c("3*mm20","3*mm30","3*mm40","3*mm21","3*mm22","3*mm31","","","","","","","3*mm10","3*mm20","3*mm30","3*mm11","3*mm12","3*mm21","","","","","","","","","","","","","","","","","","")
,c("4*mm30","4*mm40","4*mm50","4*mm31","4*mm32","4*mm41","","","","","","","6*mm20","6*mm30","6*mm40","6*mm21","6*mm22","6*mm31","","","","","","","","","","","","","","","","","","")
,c("","","","","","","1*mm00","1*mm10","1*mm20","1*mm01","1*mm02","1*mm11","","","","","","","","","","","","","","","","","","","","","","","","")
,c("","","","","","","2*mm01","2*mm11","2*mm21","2*mm02","2*mm03","2*mm12","","","","","","","","","","","","","","","","","","","1*mm00","1*mm10","1*mm20","1*mm01","1*mm02","1*mm11")
,c("","","","","","","3*mm02","3*mm12","3*mm22","3*mm03","3*mm04","3*mm13","","","","","","","","","","","","","","","","","","","3*mm01","3*mm11","3*mm21","3*mm02","3*mm03","3*mm12")
,c("","","","","","","4*mm03","4*mm13","4*mm23","4*mm04","4*mm05","4*mm14","","","","","","","","","","","","","","","","","","","6*mm02","6*mm12","6*mm22","6*mm03","6*mm04","6*mm13")
,c("1*mm01","1*mm11","1*mm21","1*mm02","1*mm03","1*mm12","1*mm10","1*mm20","1*mm30","1*mm11","1*mm12","1*mm21","","","","","","","1*mm00","1*mm10","1*mm20","1*mm01","1*mm02","1*mm11","1*mm00","1*mm10","1*mm20","1*mm01","1*mm02","1*mm11","","","","","","")
,c("1*mm02","1*mm12","1*mm22","1*mm03","1*mm04","1*mm13","2*mm11","2*mm21","2*mm31","2*mm12","2*mm13","2*mm22","","","","","","","2*mm01","2*mm11","2*mm21","2*mm02","2*mm03","2*mm12","2*mm01","2*mm11","2*mm21","2*mm02","2*mm03","2*mm12","1*mm10","1*mm20","1*mm30","1*mm11","1*mm12","1*mm21")
,c("2*mm11","2*mm21","2*mm31","2*mm12","2*mm13","2*mm22","1*mm20","1*mm30","1*mm40","1*mm21","1*mm22","1*mm31","1*mm01","1*mm11","1*mm21","1*mm02","1*mm03","1*mm12","2*mm10","2*mm20","2*mm30","2*mm11","2*mm12","2*mm21","2*mm10","2*mm20","2*mm30","2*mm11","2*mm12","2*mm21","","","","","","")
,c("2*mm12","2*mm22","2*mm32","2*mm13","2*mm14","2*mm23","2*mm21","2*mm31","2*mm41","2*mm22","2*mm23","2*mm32","1*mm02","1*mm12","1*mm22","1*mm03","1*mm04","1*mm13","4*mm11","4*mm21","4*mm31","4*mm12","4*mm13","4*mm22","4*mm11","4*mm21","4*mm31","4*mm12","4*mm13","4*mm22","1*mm20","1*mm30","1*mm40","1*mm21","1*mm22","1*mm31")
,c("1*mm03","1*mm13","1*mm23","1*mm04","1*mm05","1*mm14","3*mm12","3*mm22","3*mm32","3*mm13","3*mm14","3*mm23","","","","","","","3*mm02","3*mm12","3*mm22","3*mm03","3*mm04","3*mm13","3*mm02","3*mm12","3*mm22","3*mm03","3*mm04","3*mm13","3*mm11","3*mm21","3*mm31","3*mm12","3*mm13","3*mm22")
,c("3*mm21","3*mm31","3*mm41","3*mm22","3*mm23","3*mm32","1*mm30","1*mm40","1*mm50","1*mm31","1*mm32","1*mm41","3*mm11","3*mm21","3*mm31","3*mm12","3*mm13","3*mm22","3*mm20","3*mm30","3*mm40","3*mm21","3*mm22","3*mm31","3*mm20","3*mm30","3*mm40","3*mm21","3*mm22","3*mm31","","","","","","")
)

 namess = c('a00','a10','a20','a01','a02','a11',
            'b00','b10','b20','b01','b02','b11',
            'c00','c10','c20','c01','c02','c11',
            'd00','d10','d20','d01','d02','d11',
            'e00','e10','e20','e01','e02','e11',
            'f00','f10','f20','f01','f02','f11')
 objlist=objects(pos=1)
 checknames = rep(0,length(namess))
 for(i in 1:length(namess))
 {
    if(sum(objlist==namess[i])==1){checknames[i] = 1}
 }
 whichnames = which(checknames==1)
 func.list.timehomo=rep(0,length(namess))
 for(i in whichnames)
 {
     # which expressions vary over time
     result=eval(body(namess[i]))
     func.list.timehomo[i]=2-(sum(diff(result)==0)==(length(result)-1))
 }
 if(any(func.list.timehomo==2)){homo=F}
 dims= rep('',14*2)
 for(i in 1:14)
 {
   for(j in whichnames)
   {
     if(MAT[i,j]!='')
     {
      dims[i] =paste0(dims[i],'+(',body(namess[j])[2],')',c('*','%')[func.list.timehomo[j]],'(',MAT[i,j],')')
      dims[i+14] =paste0(dims[i+14],'+(',body(namess[j])[2],')',c('*','%')[func.list.timehomo[j]],'(',MAT2[i,j],')')

     }
   }
 }



 #namess2 = c('Lam00','Lamy0','Nmu11','Nsig11','Nmu21','Nsig21','Nmu12','Nsig12','Nmu22','Nsig22','Lam10','Lamx2','Lamy1','Lamy2','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')
 namess2 = c('Lam00','Lam10','Lam01','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')

 checknames2 = rep(0,length(namess2))
 for(i in 1:length(namess2))
 {
    if(sum(objlist==namess2[i])==1){checknames2[i] = 1}
 }
 whichnames2 = which(checknames2==1)

 func.list.timehomo2=rep(0,length(namess2))

 for(i in whichnames2)
 {
     # which expressions vary over time
     result=eval(body(namess2[i]))
     func.list.timehomo2[i]=2-(sum(diff(result)==0)==(length(result)-1))
 }
 if(any(func.list.timehomo2==2)){homo=F}
 #print(func.list.timehomo)
 #print(func.list.timehomo2)


 if(checknames2[1]==0){Lam00=function(t){0}}
 if(checknames2[2]==0){Lam10=function(t){0}}
 if(checknames2[3]==0){Lam01=function(t){0}}
 #if(checknames2[4]==0){Nsig11=function(t){0};Nmu11=function(t){0}}
 #if(checknames2[5]==0){Nmu21=function(t){0}}
 #if(checknames2[6]==0){Nsig21=function(t){0};Nmu21=function(t){0}}
 #if(checknames2[7]==0){Nmu12=function(t){0}}
 #if(checknames2[8]==0){Nsig12=function(t){0};Nmu12=function(t){0}}
 #if(checknames2[9]==0){Nmu22=function(t){0}}
 #if(checknames2[10]==0){Nsig22=function(t){0};Nmu22=function(t){0}}
 #if(checknames2[11]==0){Lam10=function(t){0}}
 #if(checknames2[12]==0){Lamx2=function(t){0}}
 #if(checknames2[13]==0){Lamy1=function(t){0}}
 #if(checknames2[14]==0){Lamy2=function(t){0}}
 if(checknames2[4]==0){Jmu1=function(t){0}}
 if(checknames2[5]==0){Jmu2=function(t){0}}
 if(checknames2[6]==0){Jsig11=function(t){0}}
 if(checknames2[7]==0){Jsig12=function(t){0}}
 if(checknames2[8]==0){Jsig22=function(t){0}}

 if(checknames2[1]==1)
 {
  for(i in 1:14)
  {
     dims[i] = paste0(dims[i],'+(',body(namess2[1])[2],')',c('*','%')[func.list.timehomo2[1]],'(',Mstar1[i],')')
  }
 }

 if(checknames2[2]==1)
 {
  for(i in 1:14)
  {
     dims[i] = paste0(dims[i],'+(',body(namess2[2])[2],')',c('*','%')[func.list.timehomo2[2]],'(',Mstar2[i],')')
  }
 }
  if(checknames2[3]==1)
 {
  for(i in 1:14)
  {
     dims[i] = paste0(dims[i],'+(',body(namess2[3])[2],')',c('*','%')[func.list.timehomo2[3]],'(',Mstar3[i],')')
  }
 }
 # if(checknames2[11]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[11])[2],')',c('*','%')[func.list.timehomo2[11]],'(',Mstar1x[i],')')
 # }
 #}

 #  if(checknames2[12]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[12])[2],')',c('*','%')[func.list.timehomo2[12]],'(',Mstar2x[i],')')
 # }
 #}
 #  if(checknames2[13]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[13])[2],')',c('*','%')[func.list.timehomo2[13]],'(',Mstar1y[i],')')
 # }
 #}

 #  if(checknames2[14]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[14])[2],')',c('*','%')[func.list.timehomo2[14]],'(',Mstar2y[i],')')
 # }
 #}
 #if((Jdist=='Normal'))
 #{
 #prem =
 #'
 #   double mm1 = mu11                                         ;
 #   double mm2 = pow(mu11,2)+ pow(sig11,2)                          ;
 #   double mm3 = pow(mu11,3)+ 3* pow(mu11,1)*pow(sig11,2)              ;
 #   double mm4 = pow(mu11,4)+ 6* pow(mu11,2)*pow(sig11,2)+3*pow(sig11,4)  ;
 #
 #   double oo1 = mu21                                         ;
 #   double oo2 = pow(mu21,2)+ pow(sig21,2)                          ;
 #   double oo3 = pow(mu21,3)+ 3* pow(mu21,1)*pow(sig21,2)              ;
 #   double oo4 = pow(mu21,4)+ 6* pow(mu21,2)*pow(sig21,2)+3*pow(sig21,4)  ;
 #
 #   double nn1 = mu12                                         ;
 #   double nn2 = pow(mu12,2)+ pow(sig12,2)                          ;
 #  double nn3 = pow(mu12,3)+ 3* pow(mu12,1)*pow(sig12,2)              ;
 #   double nn4 = pow(mu12,4)+ 6* pow(mu12,2)*pow(sig12,2)+3*pow(sig12,4)  ;

 #   double pp1 = mu22                                         ;
 #   double pp2 = pow(mu22,2)+ pow(sig22,2)                          ;
 #   double pp3 = pow(mu22,3)+ 3* pow(mu22,1)*pow(sig22,2)              ;
 #   double pp4 = pow(mu22,4)+ 6* pow(mu22,2)*pow(sig22,2)+3*pow(sig22,4)  ;
 #'
 #}

 if(Jdist=='MVNormal')
 {
  #print('MVNormal')
  prem =
 '
    double mv10 = mu1                                                          ;
    double mv20 = pow(mu1,2)+ sig11                                            ;
    double mv30 = pow(mu1,3)+ 3* pow(mu1,1)*sig11                              ;
    double mv40 = pow(mu1,4)+ 6* pow(mu1,2)*sig11+3*pow(sig11,2)               ;
    double mv01 = mu2                                                          ;
    double mv02 = pow(mu2,2)+ sig22                                            ;
    double mv03 = pow(mu2,3)+ 3* pow(mu2,1)*sig22                              ;
    double mv04 = pow(mu2,4)+ 6* pow(mu2,2)*sig22+3*pow(sig22,2)               ;
    double mv11 = mu1*mu2+sig12                                                ;
    double mv12 = mu1*pow(mu2,2)+2*mu2*sig12+mu1*sig22                         ;
    double mv21 = pow(mu1,2)*mu2+2*mu1*sig12+mu2*sig11                         ;
    double mv22 = pow(mu1,2)*pow(mu2,2)+pow(mu2,2)*sig11+pow(mu1,2)*sig22+4*mu1*mu2*sig12 +sig11*sig22+2*sig12*sig12;
    double mv13 = mu1*pow(mu2,3)+3*pow(mu2,2)*sig12+3*mu1*mu2*sig22 +3*sig12*sig22;
    double mv31 = mu2*pow(mu1,3)+3*pow(mu1,2)*sig12+3*mu1*mu2*sig11 +3*sig12*sig11;
 '
   }

     premprem =rep('',4)
     #'Lam00','Lamy0','Nmu11','Nsig11','Nmu21','Nsig21','Nmu12','Nsig12','Nmu22','Nsig22','Lam10','Lamx2','Lamy1','Lamy2'
     #if(checknames2[3]==1){premprem[1]= paste0('double mu11=',body('Nmu11')[2],';','double sig11=',body('Nsig11')[2],';')}
     #if(checknames2[5]==1){premprem[3]= paste0('double mu21=',body('Nmu21')[2],';','double sig21=',body('Nsig21')[2],';')}
     #if(checknames2[7]==1){premprem[2]= paste0('double mu12=',body('Nmu12')[2],';','double sig12=',body('Nsig12')[2],';')}
     #if(checknames2[9]==1){premprem[4]= paste0('double mu22=',body('Nmu22')[2],';','double sig22=',body('Nsig22')[2],';')}
     #namess2 = c('Lam00','Lam10','Lam01','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')
     #if(checknames2[3]==0){premprem[1]= paste0('double mu11=',0,';','double sig11=',0,';')}
     #if(checknames2[5]==0){premprem[3]= paste0('double mu21=',0,';','double sig21=',0,';')}
     #if(checknames2[7]==0){premprem[2]= paste0('double mu12=',0,';','double sig12=',0,';')}
     #if(checknames2[9]==0){premprem[4]= paste0('double mu22=',0,';','double sig22=',0,';')}
     if(Jdist=='MVNormal')
     {
        #print('MVNormal')
        #namess2 = c('Lam00','Lam10','Lam01','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')
        premprem[1]= paste0('double mu1=',body('Jmu1')[2],';\n','double mu2=',body('Jmu2')[2],';\n','double sig11=',body('Jsig11')[2],';\n','double sig12=',body('Jsig12')[2],';\n','double sig22=',body('Jsig22')[2],';\n')
        premprem[2:4]=''
     }

    prm=''
   for(i in 1:4)
   {
     prm = paste0(prm,premprem[i],'\n ')
   }

   b = rep('(+0*a.col(0))',6)
   if(checknames2[1]==1){b[1] = paste0('(',body('Lam00')[2],')',c('*','%')[func.list.timehomo2[1 ]],'m00')}else{b[1] = paste0('(+0*a.col(0))')}
   if(checknames2[2]==1){b[2] = paste0('(',body('Lam10')[2],')',c('*','%')[func.list.timehomo2[2 ]],'m10')}else{b[2] = paste0('(+0*a.col(0))')}
   if(checknames2[3]==1){b[3] = paste0('(',body('Lam01')[2],')',c('*','%')[func.list.timehomo2[3 ]],'m01')}else{b[3] = paste0('(+0*a.col(0))')}
   #if(checknames2[12]==1){b[4] = paste0('(',body('Lamx2')[2],')',c('*','%')[func.list.timehomo2[12]],'mm10')}else{b[4] = paste0('(+0*a.col(0))')}
   #if(checknames2[13]==1){b[5] = paste0('(',body('Lamy1')[2],')',c('*','%')[func.list.timehomo2[13]],'mm01')}else{b[5] = paste0('(+0*a.col(0))')}
   #if(checknames2[14]==1){b[6] = paste0('(',body('Lamy2')[2],')',c('*','%')[func.list.timehomo2[14]],'mm01')}else{b[6] = paste0('(+0*a.col(0))')}

   for(i in 1:28)
   {
      dims[i]= paste0( 'atemp.col(',i-1,')=',dims[i],';')
   }
      for(i in 1:6)
   {
      b[i]= paste0( 'atemp.col(',28+i-1,')=',b[i],';')
   }

   pr="
   vec m00 = (1+0*a.col(0));
   vec m10 = a.col(0);
   vec m20 = a.col(1);
   vec m30 = a.col(2);
   vec m40 = a.col(3);
   vec m01 = a.col(4);
   vec m02 = a.col(5);
   vec m03 = a.col(6);
   vec m04 = a.col(7);
   vec m11 = a.col(8);
   vec m12 = a.col(9);
   vec m21 = a.col(10);
   vec m22 = a.col(11);
   vec m13 = a.col(12);
   vec m31 = a.col(13);

   vec mm00 = (1+0*a.col(0));
   vec mm10 = a.col(14);
   vec mm20 = a.col(15);
   vec mm30 = a.col(16);
   vec mm40 = a.col(17);
   vec mm01 = a.col(18);
   vec mm02 = a.col(19);
   vec mm03 = a.col(20);
   vec mm04 = a.col(21);
   vec mm11 = a.col(22);
   vec mm12 = a.col(23);
   vec mm21 = a.col(24);
   vec mm22 = a.col(25);
   vec mm13 = a.col(26);
   vec mm31 = a.col(27);
   "

   pr2="
   vec m00 = (1+0*x0.col(0));
   vec m10 = x0.col(0);
   vec m20 = x0.col(1);
   vec m30 = x0.col(2);
   vec m40 = x0.col(3);
   vec m01 = x0.col(4);
   vec m02 = x0.col(5);
   vec m03 = x0.col(6);
   vec m04 = x0.col(7);
   vec m11 = x0.col(8);
   vec m12 = x0.col(9);
   vec m21 = x0.col(10);
   vec m22 = x0.col(11);
   vec m13 = x0.col(12);
   vec m31 = x0.col(13);

   vec mm00 = (1+0*x0.col(0));
   vec mm01 = x0.col(14);
   vec mm02 = x0.col(15);
   vec mm03 = x0.col(16);
   vec mm04 = x0.col(17);
   vec mm10 = x0.col(18);
   vec mm20 = x0.col(19);
   vec mm30 = x0.col(20);
   vec mm40 = x0.col(21);
   vec mm11 = x0.col(22);
   vec mm12 = x0.col(23);
   vec mm21 = x0.col(24);
   vec mm22 = x0.col(25);
   vec mm13 = x0.col(26);
   vec mm31 = x0.col(27);
   "

   odekernel=paste0(pr,prm,prem,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],
                   '\n',dims[5],'\n',dims[6],'\n',dims[7],'\n',dims[8],'\n',dims[9],'\n',dims[10],
                   '\n',dims[11],'\n',dims[12],'\n',dims[13],'\n',dims[14],'\n',dims[1+14],'\n',dims[2+14],'\n',dims[3+14],'\n',dims[4+14],
                   '\n',dims[5+14],'\n',dims[6+14],'\n',dims[7+14],'\n',dims[8+14],'\n',dims[9+14],'\n',dims[10+14],
                   '\n',dims[11+14],'\n',dims[12+14],'\n',dims[13+14],'\n',dims[14+14],'\n',b[1],'\n',b[2],'\n',b[3],'\n',b[4],'\n',b[5],'\n',b[6])


   # write.table(odekernel,'Check_ODE_dims.txt')
   # WRIGHT AND SOURCE -----------------------------------------------------------
     txt.full=paste(txtA,odekernel,txtB,txtC)
     type.sol ="                  GENERALIZED QUADRATIC DIFFUSON"
}

   #library(Rcpp)
   #library(RcppArmadillo)
   if(wrt)
   {
     write(txt.full,'BiJGQD.mcmc.cpp')
   }
   stre="Compiling C++ code. Please wait."
   cat(stre, " \r")
   flush.console()
   sourceCpp(code=txt.full)
   cat('                                     ','\r')

   seq1=c(1,2,3,4,0,0,0,0,1,1,2,2,1,3)
   seq2=c(0,0,0,0,1,2,3,4,1,2,1,2,3,1)


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
   indnames2 = rep(0,5)
   for(j in whichnames2)
   {
       if(sum(obs==namess2[j]))
       {
         indnames2[j]=TRUE
         namess2[j]=paste0(namess2[j],' : ',trim(body(namess2[j])[2]))
       }
   }
   prior.list = trim(prior.list)
   namess4=matrix(namess4,length(namess4),1)
   buffer0=c('================================================================')
   buffer1=c('----------------------------------------------------------------')
   buffer2=c('................................................................')
   buffer3=c('...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ')
   buffer4=c('_____________________ Drift Coefficients _______________________')
   buffer5=c('___________________ Diffusion Coefficients _____________________')
   buffer6=c('_____________________ Prior Distributions ______________________')
   buffer7=c('_______________________ Model/Chain Info _______________________')
   buffer8=c('......................... Intensity ............................')
   buffer9=c('........................... Jumps ..............................')
   buffer10=c('_______________________ Jump Components ________________________')
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
         namess4[31:36][which(indnames[31:36]==T)],
         buffer10,
         buffer8,
         namess2[c(1:3)][which(indnames2[1:3]==T)],
         buffer9,
         namess2[c(4:8)][which(indnames2[4:8]==T)],
         buffer6,'',prior.list)
   Info=data.frame(matrix(Info,length(Info),1))
   colnames(Info)=''
   if(print.output)
   {
     print(Info,row.names = FALSE,right=F)
   }
    ############################################################################
    ############################################################################
    ############################################################################

    tme=Sys.time()
    par.matrix=matrix(0,length(theta),updates)
    probs=rep(0,updates)
    ll=rep(0,updates)
    errors =ll
    acc=ll
    kk=0
    par.matrix[,1]=theta
    prop.matrix =par.matrix
    retries = 0
    max.retries = 0
    retry.count   = 0
    retry.indexes = c()
    success = TRUE
    strts =matrix(0,nnn-1,2)
    strts[,1]=rtf[1]
    strts[,2]=rtf[2]
    rs=solver(X1[-nnn],X2[-nnn],X1[-1],X2[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],strts,tro,secmom,seq1,seq2)
    if(is.na(sum(rs$like[,1])))
    {
          retry.count =1
          while(is.na(sum(rs$like[,1]))&&(retry.count<=10))
          {
            theta.start=theta+rnorm(length(theta),sd=sds)
            rs = solver(X1[-nnn],X2[-nnn],X1[-1],X2[-1],c(0,theta.start),mesh,delt,nnn-1,T.seq[-nnn],strts,tro,secmom,seq1,seq2)
            if(is.na(sum(rs$like[,1]))){retry.count=retry.count+1}
          }
    }


    lold=sum(rs$like[,1])
    # print(lold)
    errors[1] = rs$max
    if(is.na(lold)){print('Fail: Could not evaluate likelihood at initial perameters.');failed.chain=T;}
    ll[1]=lold
    muvec = theta
    covvec = diag(1/(adapt^2/length(theta))*(sds)^2)
    if(adapt==0)
    {
    pb <- txtProgressBar(1,updates,1,style = 1,width = 65)
    failed.chain=F
    i=2
    while(i<=updates)
    {
        theta.temp=theta
        theta=theta+rnorm(length(theta),sd=sds)
        prop.matrix[,i] = theta
        tempp=solver(X1[-nnn],X2[-nnn],X1[-1],X2[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],strts,tro,secmom,seq1,seq2)
        strts=tempp$like[,2:3]*recycle+strts*(1-recycle)
        errors[i] = tempp$max

        lnew=sum(tempp$like[,1])

        rat=min(exp(lnew-lold)*pp(theta)/pp(theta.temp),1)
        if(is.na(rat))
        {
          retry.count =1
          retries=retries+1
          retry.indexes[retries] = i
          max.retries=max.retries+1


          while(is.na(rat)&&(retry.count<=10))
          {
            i = max(i-10,2)
            theta = par.matrix[,i]
            theta=theta+rnorm(length(theta),sd=sds)
            prop.matrix[,i] = theta
            tempp=solver(X1[-nnn],X2[-nnn],X1[-1],X2[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],strts,tro,secmom,seq1,seq2)
            strts=tempp$like[,2:3]*recycle+strts*(1-recycle)
            errors[i] = tempp$max
            lnew=sum(tempp$like[,1])

            rat=min(exp(lnew-lold)*pp(theta)/pp(theta.temp),1)
            if(is.na(rat)){retry.count=retry.count+1}
            if(retry.count>10){print('Fail: Local retry fail!');failed.chain=T;break;break;}
          }
        }
        u=runif(1)
        is.true =(rat>u)
        is.false=!is.true
        theta=theta*is.true+theta.temp*is.false
        lold=lnew*is.true +lold*is.false
        par.matrix[,i]=theta
        ll[i]=lold
        probs[i] = mean(rs$probs)
        kk=kk+is.true
        acc[i]=kk/i
        if(max.retries>5000){print('Fail: Failed evaluation limit exceeded!');failed.chain=T;break;}
        if(any(is.na(theta))){print('Fail: Samples were NA! ');failed.chain=T;break;}
        setTxtProgressBar(pb, i)
        i=i+1
    }
    close(pb)
    }
    if(adapt!=0)
    {
      pb <- txtProgressBar(1,updates,1,style = 1,width = 65)
    failed.chain=F
    for(i in 2:updates)
    {
        theta.temp=theta
        if(i>min(5000,round(burns/2)))
        {
          theta=theta+(1-adapt)*rnorm(length(theta),sd=sqrt(2.38^2/length(theta)*diag(covvec)))+adapt*rnorm(length(theta),sd=sqrt(0.1^2/length(theta)*diag(covvec)))
        }else
        {
          theta=theta+rnorm(length(theta),sd=sds)
        }
        prop.matrix[,i] = theta
        tempp=solver(X1[-nnn],X2[-nnn],X1[-1],X2[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],strts,tro,secmom)
        strts=tempp$like[,2:3]
        errors[i] = tempp$max
        lnew=sum(tempp$like[-excl,1])
        rat=min(exp(lnew-lold)*pp(theta)/pp(theta.temp),1)
        u=runif(1)
        is.true =(rat>u)
        is.false=!is.true
        theta=theta*is.true+theta.temp*is.false
        lold=lnew*is.true +lold*is.false
        par.matrix[,i]=theta
        ll[i]=lold
        kk=kk+is.true
        acc[i]=kk/i
        muvec=muvec +1/(i)*(theta-muvec)
        covvec = covvec+1/(i)*((theta-muvec)%*%t(theta-muvec)-covvec)

        if(any(is.na(theta))){print('Fail: NaN thetas observed.');failed.chain=T;plot.chain=F;break;}
        setTxtProgressBar(pb, i)

    }
    close(pb)
    }
    tme.eval = function(start_time)
    {
      start_time = as.POSIXct(start_time)
      dt = difftime(Sys.time(), start_time, units="secs")
      format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
     }
    tme=tme.eval(tme)


    theta =apply(par.matrix[,-c(1:burns)],1,mean)
    meanD=mean(-2*ll[-c(1:burns)])
    rs=solver(X1[-nnn],X2[-nnn],X1[-1],X2[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],strts,tro,secmom,seq1,seq2)
    pd=meanD-(-2*sum(rs$like[-excl,1]))
    DIC=meanD+pd
    actual.p=length(theta)

    model.inf=list(elapsed.time=tme,time.homogeneous=c('Yes','No')[2-homo],p=actual.p,DIC=DIC,pd=pd,N=length(X[,1])-length(excl)+1,Tag=Tag,burns=burns)
    Info2=c( buffer7,
             paste0("Chain Updates       : ",updates),
             paste0("Burned Updates      : ",burns),
             paste0("Time Homogeneous    : ",c('Yes','No')[2-homo]),
             paste0("Data Resolution     : ",c(paste0('Homogeneous: dt=',round(max(diff(T.seq)[-excl]),4)),paste0('Variable: min(dt)=',round(min(diff(T.seq)[-excl]),4),', max(dt)=',round(max(diff(T.seq)[-excl]),4)))[2-homo.res]),
             paste0("# Removed Transits. : ",c("None",length(excl))[2-is.null(exclude)]),
             paste0("Density approx.     : ",sol.state),
             paste0('Elapsed time        : ',tme),
             buffer3,
             paste0("dim(theta)          : ",round(actual.p,3)),
             paste0("DIC                 : ",round(DIC,3)),
             paste0("pd (eff. dim(theta)): ",round(pd,3)),
             buffer1)
    Info2=data.frame(matrix(Info2,length(Info2),1))
    colnames(Info2)=''
    if(print.output)
    {
      print(Info2,row.names = FALSE,right=F)
    }
    if(plot.chain)
    {
      nper=length(theta)
      if(nper==1){par(mfrow=c(1,2))}
      if(nper==2){par(mfrow=c(2,2))}
      if(nper==3){par(mfrow=c(2,2))}
      if(nper>3)
      {
        d1=1:((nper)+1)
        d2=d1
        O=outer(d1,d2)
        test=O-((nper)+1)
        test[test<0]=100
        test=test[1:4,1:4]
        test
        wh=which(test==min(test))

        d1=d1[col(test)[wh[1]]]
        d2=d2[row(test)[wh[1]]]
        par(mfrow=c(d1,d2))
      }
      cols=rainbow_hcl(nper, start = 10, end = 275,c=100,l=70)
      ylabs=paste0('theta[',1:nper,']')
      for(i in 1:nper)
      {
          plot(prop.matrix[i,],col='gray90',type='s',main=ylabs[i],xlab='Iteration',ylab='')
          lines(par.matrix[i,],col=cols[i],type='s')
          abline(v=burns,lty='dotdash')
          if(adapt!=0){abline(v=min(5000,round(burns/2)),lty='dotted',col='red')}
      }
      plot(acc,type='l',ylim=c(0,1),col='darkblue',main='Accept. Rate',xlab='Iteration',ylab='%/100')
      abline(h=seq(0,1,1/10),lty='dotted')
      abline(v=burns,lty='dotdash')
      abline(h=0.4,lty='solid',col='red',lwd=1.2)
      abline(h=0.2,lty='solid',col='red',lwd=1.2)
      lines(probs,col='purple',lty='dotted',lwd=2)

      box()
    }
    ret=list(par.matrix=t(par.matrix),acceptance.rate=acc,elapsed.time=tme,model.info=model.inf,failed.chain=failed.chain,covvec=covvec,prop.matrix=t(prop.matrix),errors=errors)
    class(ret) = 'JGQD.mcmc'
    return(ret)
  }




