GQD.mle <-function(X,time,mesh=10,theta,control=NULL,method='Nelder-Mead',Dtype='Saddle',Trunc=c(4,4),RK.order=4,P=200,alpha=0,lower=min(X)/2,upper=max(X)*2,exclude=NULL,Tag=NA,wrt=FALSE,print.output=TRUE)
{
  solver   =function(Xs, Xt, theta, N , delt , N2, tt  , P , alpha, lower , upper, tro  ){}
  rm(list =c('solver'))  
  theta = theta+runif(length(theta),0.01,0.02)*sign(theta)
     check_for_model=function()
  {
    txt=''
    namess=c('G0','G1','G2','Q0','Q1','Q2')
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
  G0=function(t){theta[1]*theta[2]}
  G1=function(t){-theta[1]}
  Q1=function(t){theta[3]*theta[3]}

  model=GQD.mle(X,time,10,theta =rep(1,3))
  --------------------------------------------------------------------------------
   '
   check=T
   }
   if((sum(func.list)>0)&&(sum(func.list[-c(1:3)])==0))
   {
   txt='
  --------------------------------------------------------------------------------
  At least one diffusion coefficient has to be defined! Try for example:
  --------------------------------------------------------------------------------
  GQD.remove()
  G0=function(t){theta[1]*theta[1]}
  model=GQD.mle(X,time,10,theta =c(0.1))
  --------------------------------------------------------------------------------
   '
     check=T
    }
    return(list(check=check,txt=txt))
  }
  check_for=check_for_model()
  if(check_for[[1]]){stop(check_for[[2]])}

  theta = theta+runif(length(theta),0.001,0.002)*sign(theta)
  pow=function(x,p)
  {
    x^p
  }
  prod=function(a,b){a*b}
  T.seq=time 
  # Warning Module
  TR.order=Trunc[1]
  DTR.order=Trunc[2]
  Dtypes =c('Saddle','Normal','Gamma','InvGamma','Beta')
  Dindex = which(Dtypes==Dtype)
  IntRange = c(lower,upper)
  IntDelta =1/100
  Xtr = min(IntRange)
     b1 = '\n==============================================================================\n'
   b2 = '==============================================================================\n'
  warn=c(
    '1.  Missing input: Argument {X} is missing.\n'
    ,'2.  Missing input: Argument {time} is missing.\n'
    ,'3.  Missing input: Argument {theta} is missing.\n'
    ,'4.  Missing input: Argument {sds} is missing.\n'
    ,'5.  Input type: Argument {X} must be of type vector!.\n'
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
    ,'18. Density: Dtype has to be one of Saddle, Normal, Gamma, InvGamma or Beta.\n'
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
     namess=c('G0','G1','G2','Q0','Q1','Q2')
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
         theta[j] = theta[j]+runif(1,-0.01,0.01)
         dresult2=eval(body(namess[i]))
         dff = abs(dresult1-dresult2)
         if(any(round(dff,6)!=0)){pers.represented[j]=pers.represented[j]+1}
       }
     }
     return(pers.represented)
   }

   check.thetas2 = function(theta)
   {
     namess=c('G0','G1','G2','Q0','Q1','Q2')
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
   if(!is.vector(X))                                             {warntrue[5]=T}
   if(!is.vector(time))                                          {warntrue[6]=T}
   # Check model parameters:
   if(check.thetas2(theta)!=0)                                   {warntrue[7]=T}
   if(!warntrue[7]){if(any(check.thetas(theta,T.seq)==0))        {warntrue[8]=T}}

   # Check input length:
   if(length(X)<10)                                              {warntrue[9]=T}
   if(length(time)<10)                                          {warntrue[10]=T}
   if(length(lower)>1)                                          {warntrue[11]=T}
   if(length(upper)>1)                                          {warntrue[12]=T}
   if(length(P)!=1)                                             {warntrue[13]=T}
   if(length(mesh)!=1)                                        {warntrue[14]=T}
   if(length(alpha)!=1)                                         {warntrue[15]=T}
   if(length(Trunc)!=2)                                         {warntrue[16]=T}
   if(length(RK.order)!=1)                                      {warntrue[17]=T}
   #if(length(N.updates)!=1)                                     {warntrue[33]=T}
   #if(length(burns)!=1)                                         {warntrue[34]=T}
   if(length(Dtype)!=1)                                         {warntrue[37]=T}


   # Check density approx parameters:
   if(sum(Dindex)==0)                                          {warntrue[18] =T}
   if(!warntrue[18])
   {
    if((Dindex==3)|(Dindex==4)){if(lower[1]<=0)                {warntrue[19] =T}}
    if(Dindex==5){if(any(X<=0)|any(X>=1))                      {warntrue[20] =T}}
   }
   if(!any(warntrue[c(11,12)])){if(upper<=lower)                {warntrue[21] =T}}
   if(!warntrue[13]){if(P<10)                                   {warntrue[22] =T}}
   if(!warntrue[16]){if(Trunc[2]>Trunc[1])                      {warntrue[23] =T}}

   #  Miscelaneous checks:
   excl=0
   if(is.null(exclude)){excl=length(T.seq)-1+200}
   if(!is.null(exclude)){excl=exclude}
   test.this =max(diff(T.seq)[-excl])/mesh
   if(test.this>0.1)                                            {warntrue[24]=T}
   if(test.this>=1)                                             {warntrue[25]=T}
   if(!warntrue[17]){if(!((RK.order==4)|(RK.order==10)))        {warntrue[26]=T}}
   if(!warntrue[14]){if(mesh<5)                               {warntrue[27]=T}}
   if(length(X)!=length(time))                                  {warntrue[28]=T}
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


   nnn=length(X)
  
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
  namess=c('G0','G1','G2','Q0','Q1','Q2')
  func.list=rep(0,6)
  obs=objects(pos=1)
  for(i in 1:6)
  {
    if(sum(obs==namess[i])){func.list[i]=1}
  }
  
  state1=(sum(func.list[c(3,5,6)]==1)==0)
  if(state1){DTR.order=2;TR.order=2;sol.state='Normally distributed diffusion.';} 
  if((state1&(Dtype!='Saddlepoint'))){TR.order=2;DTR.order=2;sol.state='2nd Ord. Truncation + Std Normal Dist.';}
  state2=!state1
  if(state2)
  {
    state2.types=c('Saddlepoint Appr.   ',
                   ' Ext. Normal Appr.    ',
                   ' Ext. Gamma Appr.     ',
                   ' Ext. Inv. Gamma Appr.')
    sol.state=paste0(TR.order,' Ord. Truncation +',DTR.order,'th Ord. ',state2.types[Dindex])
  }
  
  
  
  if(TR.order==2)
  {
    fpart=
      '#include <RcppArmadillo.h>
    #include <math.h>
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
    
    mat atemp(N2,2);'
}
if(TR.order==4)
{
  fpart=
    '#include <RcppArmadillo.h>
  #include <math.h>
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
  
  mat atemp(N2,4);'
}
if(TR.order==6)
{
  fpart=
    '
  #include <RcppArmadillo.h>
  #include <math.h>
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
  
  mat atemp(N2,6);'
}
if(TR.order==8)
{
  fpart=
    '
  #include <RcppArmadillo.h>
  #include <math.h>
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
}
if(RK.order==10)
{
  if(homo.res)
  {
    ODEpart= '
    return atemp;
  }
    
    // [[Rcpp::export]]
    vec  solver(vec Xs,vec Xt,vec theta,int N,double delt,int N2,vec tt,int P,double alpha,double lower,double upper,int tro)
{
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
    vec d=tt;
    for (int i = 1; i < N+1; i++)
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
  ODEpart= '
  return atemp;
}
  
  // [[Rcpp::export]]
  vec  solver(vec Xs,vec Xt,vec theta,int N, mat delt,int N2,vec d,int P,double alpha,double lower,double upper,int tro)
{
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
  for (int i = 1; i < N+1; i++)
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
if(RK.order==4)
{
  if(homo.res)
  {
    ODEpart= '
    return atemp;
  }
    
    // [[Rcpp::export]]
    vec  solver(vec Xs,vec Xt,vec theta,int N,double delt,int N2,vec tt,int P,double alpha,double lower,double upper,int tro)
{
    mat x0(N2,tro);
    mat fx0(N2,tro);
    mat fx1(N2,tro);
    mat fx2(N2,tro);
    mat fx3(N2,tro);
    mat x1(N2,tro);
    mat x2(N2,tro);
    mat x3(N2,tro);
    x0.fill(0);
    x0.col(0)=Xs;
    vec d=tt;
    for (int i = 1; i < N+1; i++)
{
    
    fx0 =f(x0,theta,d,N2);
    x1  =x0+delt*(0.5*fx0);
    fx1 =f(x1,theta,d+0.5*delt,N2);
    x2  =x0+delt*(0.5*fx1);
    fx2 =f(x2,theta,d+0.5*delt,N2);
    x3  =x0+delt*(fx2);
    fx3 =f(x3,theta,d+delt,N2);
    x0=x0+(0.1666667*fx0+0.3333333*fx1+0.3333333*fx2+0.1666667*fx3)*delt;
    d=d+delt;
}
    '
}
if(!homo.res)
{
  ODEpart= '
  return atemp;
}
  
  // [[Rcpp::export]]
  vec  solver(vec Xs,vec Xt,vec theta,int N, mat delt,int N2,vec d,int P,double alpha,double lower,double upper,int tro)
{
  mat x0(N2,tro);
  mat fx0(N2,tro);
  mat fx1(N2,tro);
  mat fx2(N2,tro);
  mat fx3(N2,tro);
  mat x1(N2,tro);
  mat x2(N2,tro);
  mat x3(N2,tro);
  x0.fill(0);
  x0.col(0)=Xs;
  for (int i = 1; i < N+1; i++)
{
  
  fx0 =f(x0,theta,d,N2);
  x1  =x0+delt%(0.5*fx0);
  fx1 =f(x1,theta,d+0.5*delt.col(0),N2);
  x2  =x0+delt%(0.5*fx1);
  fx2 =f(x2,theta,d+0.5*delt.col(0),N2);
  x3  =x0+delt%(fx2);
  fx3 =f(x3,theta,d+delt.col(0),N2);
  x0=x0+(0.1666667*fx0+0.3333333*fx1+0.3333333*fx2+0.1666667*fx3)%delt;
  d=d+delt.col(0);
}
  '
}
}

if(DTR.order==2)
{
  Dpart='
  return(-0.5*log(2*3.141592653589793)-0.5*log(x0.col(1))-0.5*((Xt-x0.col(0))%(Xt-x0.col(0))/x0.col(1)));
}'
}
if(DTR.order==4)
{
  Inv=
    '
  mat u(N2,4);
  u.col(0)=x0.col(0);
  u.col(1)=x0.col(1)+x0.col(0)%u.col(0);
  u.col(2)=x0.col(2)+x0.col(0)%u.col(1)+2*x0.col(1)%u.col(0);
  u.col(3)=x0.col(3)+x0.col(0)%u.col(2)+3*x0.col(1)%u.col(1)+3*x0.col(2)%u.col(0);
  
  vec det=u.col(1)%u.col(3)+u.col(0)%u.col(2)%u.col(1)+u.col(1)%u.col(0)%u.col(2)-u.col(2)%u.col(2)-u.col(1)%u.col(1)%u.col(1)-u.col(0)%u.col(0)%u.col(3);
  
  vec b11=(u.col(1)%u.col(3)-u.col(2)%u.col(2))/det;
  vec b12=(u.col(1)%u.col(2)-u.col(0)%u.col(3))/det;
  vec b13=(u.col(0)%u.col(2)-u.col(1)%u.col(1))/det;
  
  vec b21=(u.col(2)%u.col(1)-u.col(0)%u.col(3))/det;
  vec b22=(u.col(3)-u.col(1)%u.col(1))/det;
  vec b23=(u.col(1)%u.col(0)-u.col(2))/det;
  
  vec b31=(u.col(0)%u.col(2)-u.col(1)%u.col(1))/det;
  vec b32=(u.col(0)%u.col(1)-u.col(2))/det;
  vec b33=(u.col(1)-u.col(0)%u.col(0))/det;
  '
  switch(Dindex,
{
  Dpart='
  vec p=(1.0/3.0) *(3*(x0.col(3)/6.0)%x0.col(1) - pow(x0.col(2)/2.0,2))/pow(x0.col(3)/6.0,2);
  vec q=(1.0/27.0)*(27*pow(x0.col(3)/6.0,2)%(x0.col(0)-Xt) - 9*(x0.col(3)/6.0)%(x0.col(2)/2.0)%x0.col(1) + 2*pow(x0.col(2)/2.0,3))/pow(x0.col(3)/6.0,3);
  vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  vec th=-(x0.col(2)/2.0)/(3*(x0.col(3)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));
  
  vec K =x0.col(0)%th+(x0.col(1)%th%th)/2.0+(x0.col(2)%th%th%th)/6.0 +(x0.col(3)%th%th%th%th)/24.0;
  vec K1=x0.col(0)   +(x0.col(1)%th)       +(x0.col(2)%th%th)/2.0    +(x0.col(3)%th%th%th)/6.0;
  vec K2=x0.col(1)   +(x0.col(2)%th)       +(x0.col(3)%th%th)/2.0;
  vec val=-0.5*log(2*3.141592653589793*K2)+(K-th%K1);
  return(val);
}
  '
},
{
  Dpart=
    '
  vec betas1 =+b12+2*b13%u.col(0);
  vec betas2 =+b22+2*b23%u.col(0);
  vec betas3 =+b32+2*b33%u.col(0);

  vec lo =-(0.5/(u.col(0)-lower))%(sqrt(exp(2*alpha)+4*pow((u.col(0)-lower),2))-exp(alpha));
  vec up =-(0.5/(u.col(0)-upper))%(sqrt(exp(2*alpha)+4*pow((u.col(0)-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u.col(0);
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  vec K=exp(-betas1%pow(tau,1)-0.5*betas2%pow(tau,2)-0.333333333333333333*betas3%pow(tau,3))%rho%DT;
  for (int i = 1; i <= P; i++)
  {
    lo=lo+DT;
    tau   = exp(alpha)*lo/(1-pow(lo,2))+u.col(0);
    rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
    K=K+exp(-betas1%pow(tau,1)-0.5*betas2%pow(tau,2)-0.333333333333333333*betas3%pow(tau,3))%rho%DT;
  }
  return((-betas1%pow(Xt,1)-0.5*betas2%pow(Xt,2)-0.333333333333333333*betas3%pow(Xt,3))-log(K));
}
  '
},
{
  Dpart='
  vec betas1 =+b11+2*b12%u.col(0)+3*b13%u.col(1);
  vec betas2 =+b21+2*b22%u.col(0)+3*b23%u.col(1);
  vec betas3 =+b31+2*b32%u.col(0)+3*b33%u.col(1);
  
    vec lo =-(0.5/(u.col(0)-lower))%(sqrt(exp(2*alpha)+4*pow((u.col(0)-lower),2))-exp(alpha));
  vec up =-(0.5/(u.col(0)-upper))%(sqrt(exp(2*alpha)+4*pow((u.col(0)-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u.col(0);
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  vec K=exp(-betas1%log(tau)+(-betas2%tau-0.5*betas3%pow(tau,2)))%rho%DT;
  for (int i = 1; i <= P; i++)
  {
    lo=lo+DT;
    tau   = exp(alpha)*lo/(1-pow(lo,2))+u.col(0);
    rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
    K=K+exp(-betas1%log(tau)+(-betas2%tau-0.5*betas3%pow(tau,2)))%rho%DT;
  }
  return((-betas1%log(Xt)+(-betas2%Xt-0.5*betas3%pow(Xt,2)))-log(K));
}
  '
},
{
  Dpart='
  vec betas1 =+2*b11%u.col(0)+3*b12%u.col(1)+4*b13%u.col(2);
  vec betas2 =+2*b21%u.col(0)+3*b22%u.col(1)+4*b23%u.col(2);
  vec betas3 =+2*b31%u.col(0)+3*b32%u.col(1)+4*b33%u.col(2);
  
  vec lo =-(0.5/(u.col(0)-lower))%(sqrt(exp(2*alpha)+4*pow((u.col(0)-lower),2))-exp(alpha));
  vec up =-(0.5/(u.col(0)-upper))%(sqrt(exp(2*alpha)+4*pow((u.col(0)-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u.col(0);
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  vec K=exp(-betas2%log(tau)+(betas1/tau-betas3%tau))%rho%DT;
  for (int i = 1; i <= P; i++)
  {
    lo=lo+DT;
    tau   = exp(alpha)*lo/(1-pow(lo,2))+u.col(0);
    rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
    K=K+exp(-betas2%log(tau)+(betas1/tau-betas3%tau))%rho%DT;
  }
  return((-betas2%log(Xt)+(betas1/Xt-betas3%Xt))-log(K));
}
  '
},
{
  Dpart='
  
  vec betas1 =+b11%(1-2*u.col(0))+b12%(2*u.col(0)-3*u.col(1))+b13%(3*u.col(1)-4*u.col(2));
  vec betas2 =+b21%(1-2*u.col(0))+b22%(2*u.col(0)-3*u.col(1))+b23%(3*u.col(1)-4*u.col(2));
  vec betas3 =+b31%(1-2*u.col(0))+b32%(2*u.col(0)-3*u.col(1))+b33%(3*u.col(1)-4*u.col(2));
  
}'
})
if(Dindex!=1)
{
  Dpart=paste(Inv,Dpart)
}
}
if(DTR.order==6)
{
  Inv=
    '
  vec u1 = x0.col(0);
  vec u2 = x0.col(1)+x0.col(0)%u1;
  vec u3 = x0.col(2)+x0.col(0)%u2+2*x0.col(1)%u1;
  vec u4 = x0.col(3)+x0.col(0)%u3+3*x0.col(1)%u2+3*x0.col(2)%u1 ;
  vec u5 = x0.col(4)+x0.col(0)%u4+4*x0.col(1)%u3+6*x0.col(2)%u2+4*x0.col(3)%u1 ;
  vec u6 = x0.col(5)+x0.col(0)%u5+5*x0.col(1)%u4+10*x0.col(2)%u3+10*x0.col(3)%u2+5*x0.col(4)%u1    ;
  
  vec det=-u1%(u1%(u4%u6-u5%u5)-u3%(u2%u6-u3%u5)+u4%(u2%u5-u3%u4))+u2%(u1%(u3%u6-u4%u5)-u2%(u2%u6-u3%u5)+u4%(u2%u4-u3%u3))+u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)-u3%(u1%(u3%u5-u4%u4)-u2%(u2%u5-u3%u4)+u3%(u2%u4-u3%u3))+u4%(u3%u5-u4%u4);
  
  vec b11 = (u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4))/det ;
  vec b12 = (-u1%(u4%u6-u5%u5)+u2%(u3%u6-u4%u5)-u3%(u3%u5-u4%u4))/det ;
  vec b13 = (u1%(u3%u6-u4%u5)-u2%(u2%u6-u4%u4)+u3%(u2%u5-u3%u4))/det  ;
  vec b14 = (-u1%(u3%u5-u4%u4)+u2%(u2%u5-u3%u4)-u3%(u2%u4-u3%u3))/det ;
  
  vec b21 = (-u1%(u4%u6-u5%u5)+u3%(u2%u6-u3%u5)-u4%(u2%u5-u3%u4))/det;
  vec b22 = (-u2%(u2%u6-u3%u5)+u4%u6-u5%u5+u3%(u2%u5-u3%u4))/det ;
  vec b23 = (u2%(u1%u6-u3%u4)-u3%u6-u3%(u1%u5-u3%u3)+u4%u5)/det  ;
  vec b24 = (-u2%(u1%u5-u2%u4)+u3%u5-u4%u4+u3%(u1%u4-u2%u3))/det ;
  
  vec b31 = (u1%(u3%u6-u4%u5)-u2%(u2%u6-u3%u5)+u4%(u2%u4-u3%u3))/det;
  vec b32 = (u1%(u2%u6-u3%u5)-u3%u6+u4%u5-u3%(u2%u4-u3%u3))/det ;
  vec b33 = (-u1%(u1%u6-u3%u4)+u2%u6-u4%u4+u3%(u1%u4-u2%u3))/det;
  vec b34 = (u1%(u1%u5-u2%u4)-u2%u5+u3%u4-u3%(u1%u3-u2%u2))/det ;
  
  vec b41 = (-u1%(u3%u5-u4%u4)+u2%(u2%u5-u3%u4)-u3%(u2%u4-u3%u3))/det;
  vec b42 = (-u1%(u2%u5-u3%u4)+u3%u5-u4%u4+u2%(u2%u4-u3%u3))/det;
  vec b43 = (u1%(u1%u5-u3%u3)-u2%u5-u2%(u1%u4-u2%u3)+u3%u4)/det ;
  vec b44 = (-u1%(u1%u4-u2%u3)+u2%u4-u3%u3+u2%(u1%u3-u2%u2))/det ;
  '
  switch(Dindex,
{
  Dpart='
  vec p=(1.0/3.0) *(3*(x0.col(3)/6.0)%x0.col(1) - pow(x0.col(2)/2.0,2))/pow(x0.col(3)/6.0,2);
  vec q=(1.0/27.0)*(27*pow(x0.col(3)/6.0,2)%(x0.col(0)-Xt) - 9*(x0.col(3)/6.0)%(x0.col(2)/2.0)%x0.col(1) + 2*pow(x0.col(2)/2.0,3))/pow(x0.col(3)/6.0,3);
  vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  vec th=-(x0.col(2)/2.0)/(3*(x0.col(3)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));
  
  vec K =x0.col(0)%th+(x0.col(1)%th%th)/2.0+(x0.col(2)%th%th%th)/6.0 +(x0.col(3)%th%th%th%th)/24.0;
  vec K1=x0.col(0)   +(x0.col(1)%th)       +(x0.col(2)%th%th)/2.0    +(x0.col(3)%th%th%th)/6.0;
  vec K2=x0.col(1)   +(x0.col(2)%th)       +(x0.col(3)%th%th)/2.0;
  vec val=-0.5*log(2*3.141592653589793*K2)+(K-th%K1);
  return(val);
}
  '
},
{
  Dpart=
    '
  vec betas1 =+b12+2*b13%u1+3*b14%u2;
  vec betas2 =+b22+2*b23%u1+3*b24%u2;
  vec betas3 =+b32+2*b33%u1+3*b34%u2;
  vec betas4 =+b42+2*b43%u1+3*b44%u2;
  
  vec lo =-(0.5/(u1-lower))%(sqrt(exp(2*alpha)+4*pow((u1-lower),2))-exp(alpha));
  vec up =-(0.5/(u1-upper))%(sqrt(exp(2*alpha)+4*pow((u1-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  vec K=exp(-betas1%pow(tau,1)-0.5*betas2%pow(tau,2)-0.333333333333333333*betas3%pow(tau,3)-0.25*betas4%pow(tau,4))%rho%DT;
  for (int i = 1; i <= P; i++)
{
  lo=lo+DT;
  tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  K=K+exp(-betas1%pow(tau,1)-0.5*betas2%pow(tau,2)-0.333333333333333333*betas3%pow(tau,3)-0.25*betas4%pow(tau,4))%rho%DT;
}
  return((-betas1%pow(Xt,1)-0.5*betas2%pow(Xt,2)-0.333333333333333333*betas3%pow(Xt,3)-0.25*betas4%pow(Xt,4))-log(K));
}
  '
},
{
  Dpart='
  
  vec betas1 =+b11+2*b12%u1+3*b13%u2+4*b14%u3;
  vec betas2 =+b21+2*b22%u1+3*b23%u2+4*b24%u3;
  vec betas3 =+b31+2*b32%u1+3*b33%u2+4*b34%u3;
  vec betas4 =+b41+2*b42%u1+3*b43%u2+4*b44%u3;
  
  vec lo =-(0.5/(u1-lower))%(sqrt(exp(2*alpha)+4*pow((u1-lower),2))-exp(alpha));
  vec up =-(0.5/(u1-upper))%(sqrt(exp(2*alpha)+4*pow((u1-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  vec K=exp(-betas1%log(tau)+(-betas2%tau-0.5*betas3%pow(tau,2)-0.33333333333333333333333*betas4%pow(tau,3)))%rho%DT;
  for (int i = 1; i <= P; i++)
{
  lo=lo+DT;
  tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  K=K+exp(-betas1%log(tau)+(-betas2%tau-0.5*betas3%pow(tau,2)-0.33333333333333333333333*betas4%pow(tau,3)))%rho%DT;
}
  return((-betas1%log(Xt)+(-betas2%Xt-0.5*betas3%pow(Xt,2)-0.33333333333333333333333*betas4%pow(Xt,3)))-log(K));
}
  '
},
{
  Dpart='
  
  vec betas1 =+2*b11%u1+3*b12%u2+4*b13%u3+5*b14%u4;
  vec betas2 =+2*b21%u1+3*b22%u2+4*b23%u3+5*b24%u4;
  vec betas3 =+2*b31%u1+3*b32%u2+4*b33%u3+5*b34%u4;
  vec betas4 =+2*b41%u1+3*b42%u2+4*b43%u3+5*b44%u4;
  
  vec lo =-(0.5/(u1-lower))%(sqrt(exp(2*alpha)+4*pow((u1-lower),2))-exp(alpha));
  vec up =-(0.5/(u1-upper))%(sqrt(exp(2*alpha)+4*pow((u1-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  
  vec K=exp(-betas2%log(tau)+(betas1/tau-betas3%tau-0.5*betas4%pow(tau,2)))%rho%DT;
  for (int i = 1; i <= P; i++)
{
  lo=lo+DT;
  tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  K=K+exp(-betas2%log(tau)+(betas1/tau-betas3%tau-0.5*betas4%pow(tau,2)))%rho%DT;
}
  return((-betas2%log(Xt)+(betas1/Xt-betas3%Xt-0.5*betas4%pow(Xt,2)))-log(K));
}
  '
},
{
  Dpart='
  vec betas1 =+b11%(1-2*u1)+b12%(2*u1-3*u2)+b13%(3*u2-4*u3)+b14%(4*u3-5*u4);
  vec betas2 =+b21%(1-2*u1)+b22%(2*u1-3*u2)+b23%(3*u2-4*u3)+b24%(4*u3-5*u4);
  vec betas3 =+b31%(1-2*u1)+b32%(2*u1-3*u2)+b33%(3*u2-4*u3)+b34%(4*u3-5*u4);
  vec betas4 =+b41%(1-2*u1)+b42%(2*u1-3*u2)+b43%(3*u2-4*u3)+b44%(4*u3-5*u4);
  
  vec K=exp(-betas1*log(Xtr)+(betas1+betas2+betas3+betas4)*log(1-Xtr)+(betas3+betas4)*Xtr+0.5*betas4*pow(Xtr,2));
  for (int i = 1; i <= lim; i++)
{
  Xtr=Xtr+delt2;
  K=K+exp(-betas1*log(Xtr)+(betas1+betas2+betas3+betas4)*log(1-Xtr)+(betas3+betas4)*Xtr+0.5*betas4*pow(Xtr,2))*delt2;
}
  return((-betas1%log(Xt)+(betas1+betas2+betas3+betas4)%log(1-Xt)+(betas3+betas4)%Xt+0.5*betas4%pow(Xt,2))-log(K));
}
  '
})
if(Dindex!=1)
{
  Dpart=paste(Inv,Dpart)
}
}
if(DTR.order==8)
{
  Inv='
  vec u1=                                                                                  x0.col(0);
  vec u2=                                                                       x0.col(1)+1*x0.col(0)%u1;
  vec u3=                                                            x0.col(2)+1*x0.col(0)%u2+2*x0.col(1)%u1;
  vec u4=                                                 x0.col(3)+1*x0.col(0)%u3+3*x0.col(1)%u2+3*x0.col(2)%u1;
  vec u5=                                      x0.col(4)+1*x0.col(0)%u4+4*x0.col(1)%u3+6*x0.col(2)%u2+4*x0.col(3)%u1;
  vec u6=                         x0.col(5)+1*x0.col(0)%u5+5*x0.col(1)%u4+10*x0.col(2)%u3+10*x0.col(3)%u2+5*x0.col(4)%u1;
  vec u7=             x0.col(6)+1*x0.col(0)%u6+6*x0.col(1)%u5+15*x0.col(2)%u4+20*x0.col(3)%u3+15*x0.col(4)%u2+6*x0.col(5)%u1;
  vec u8= x0.col(7)+1*x0.col(0)%u7+7*x0.col(1)%u6+21*x0.col(2)%u5+35*x0.col(3)%u4+35*x0.col(4)%u3+21*x0.col(5)%u2+7*x0.col(6)%u1;
  
  vec det=-u1%(u1%(u4%(u6%u8-u7%u7)-u5%(u5%u8-u6%u7)+u6%(u5%u7-u6%u6))-u3%(u2%(u6%u8-u7%u7)-u5%(u3%u8-u4%u7)+u6%(u3%u7-u4%u6))+u4%
  (u2%(u5%u8-u6%u7)-u4%(u3%u8-u4%u7)+u6%(u3%u6-u4%u5))-u5%(u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)+u5%(u3%u6-u4%u5)))+u2%(u1%
  (u3%(u6%u8-u7%u7)-u5%(u4%u8-u5%u7)+u6%(u4%u7-u5%u6))-u2%(u2%(u6%u8-u7%u7)-u5%(u3%u8-u4%u7)+u6%(u3%u7-u4%u6))+u4%
  (u2%(u4%u8-u5%u7)-u3%(u3%u8-u4%u7)+(u3%u5-u4%u4)%u6)-u5%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u4%u6)+u5%(u3%u5-u4%u4)))-u3%(u1%
  (u3%(u5%u8-u6%u7)-u4%(u4%u8-u5%u7)+u6%(u4%u6-u5%u5))-u2%(u2%(u5%u8-u6%u7)-u4%(u3%u8-u4%u7)+u6%(u3%u6-u4%u5))+u3%
  (u2%(u4%u8-u5%u7)-u3%(u3%u8-u4%u7)+(u3%u5-u4%u4)%u6)-u5%(u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4)))+u2%
  (u4%(u6%u8-u7%u7)-u5%(u5%u8-u6%u7)+u6%(u5%u7-u6%u6))-u3%(u3%(u6%u8-u7%u7)-u5%(u4%u8-u5%u7)+u6%(u4%u7-u5%u6))+u4%
  (u3%(u5%u8-u6%u7)-u4%(u4%u8-u5%u7)+u6%(u4%u6-u5%u5))+u4%(u1%(u3%(u5%u7-u6%u6)-u4%(u4%u7-u5%u6)+u5%(u4%u6-u5%u5))-u2%
  (u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)+u5%(u3%u6-u4%u5))+u3%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u4%u6)+u5%(u3%u5-u4%u4))-u4%
  (u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4)))-u5%(u3%(u5%u7-u6%u6)-u4%(u4%u7-u5%u6)+u5%(u4%u6-u5%u5));
  
  vec b11= (u2%(u4%(u6%u8-u7%u7)-u5%(u5%u8-u6%u7)+u6%(u5%u7-u6%u6))-u3%(u3%(u6%u8-u7%u7)-u5%(u4%u8-u5%u7)+u6%(u4%u7-u5%u6))+u4%(u3%(u5%u8-u6%u7)-u4%(u4%u8-u5%u7)+u6%(u4%u6-u5%u5))-u5%(u3%(u5%u7-u6%u6)-u4%(u4%u7-u5%u6)+u5%(u4%u6-u5%u5)))/det;
  vec b12= (-u1%(u4%(u6%u8-u7%u7)-u5%(u5%u8-u6%u7)+u6%(u5%u7-u6%u6))+u2%(u3%(u6%u8-u7%u7)-u5%(u4%u8-u5%u7)+u6%(u4%u7-u5%u6))-u3%(u3%(u5%u8-u6%u7)-u4%(u4%u8-u5%u7)+u6%(u4%u6-u5%u5))+u4%(u3%(u5%u7-u6%u6)-u4%(u4%u7-u5%u6)+u5%(u4%u6-u5%u5)))/det;
  vec b13= (u1%(u3%(u6%u8-u7%u7)-u4%(u5%u8-u6%u7)+u5%(u5%u7-u6%u6))-u2%(u2%(u6%u8-u7%u7)-u4%(u4%u8-u5%u7)+u5%(u4%u7-u5%u6))+u3%(u2%(u5%u8-u6%u7)-u3%(u4%u8-u5%u7)+u5%(u4%u6-u5%u5))-u4%(u2%(u5%u7-u6%u6)-u3%(u4%u7-u5%u6)+u4%(u4%u6-u5%u5)))/det;
  vec b14= (-u1%(u3%(u5%u8-u6%u7)-u4%(u4%u8-u6%u6)+u5%(u4%u7-u5%u6))+u2%(u2%(u5%u8-u6%u7)-u4%(u3%u8-u5%u6)+u5%(u3%u7-u5%u5))-u3%(u2%(u4%u8-u6%u6)-u3%(u3%u8-u5%u6)+u5%(u3%u6-u4%u5))+u4%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u5%u5)+u4%(u3%u6-u4%u5)))/det;
  vec b15= (u1%(u3%(u5%u7-u6%u6)-u4%(u4%u7-u5%u6)+u5%(u4%u6-u5%u5))-u2%(u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)+u5%(u3%u6-u4%u5))+u3%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u4%u6)+u5%(u3%u5-u4%u4))-u4%(u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4)))/det;
  
  vec b21= (-u1%(u4%(u6%u8-u7%u7)-u5%(u5%u8-u6%u7)+u6%(u5%u7-u6%u6))+u3%(u2%(u6%u8-u7%u7)-u5%(u3%u8-u4%u7)+u6%(u3%u7-u4%u6))-u4%(u2%(u5%u8-u6%u7)-u4%(u3%u8-u4%u7)+u6%(u3%u6-u4%u5))+u5%(u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)+u5%(u3%u6-u4%u5)))/det;
  vec b22= (-u2%(u2%(u6%u8-u7%u7)-u5%(u3%u8-u4%u7)+u6%(u3%u7-u4%u6))+u3%(u2%(u5%u8-u6%u7)-u4%(u3%u8-u4%u7)+u6%(u3%u6-u4%u5))+u4%(u6%u8-u7%u7)-u5%(u5%u8-u6%u7)-u4%(u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)+u5%(u3%u6-u4%u5))+u6%(u5%u7-u6%u6))/det;
  vec b23= (u2%(u1%(u6%u8-u7%u7)-u4%(u3%u8-u4%u7)+u5%(u3%u7-u4%u6))-u3%(u1%(u5%u8-u6%u7)-u3%(u3%u8-u4%u7)+u5%(u3%u6-u4%u5))-u3%(u6%u8-u7%u7)+u4%(u5%u8-u6%u7)+u4%(u1%(u5%u7-u6%u6)-u3%(u3%u7-u4%u6)+u4%(u3%u6-u4%u5))-u5%(u5%u7-u6%u6))/det;
  vec b24= (-u2%(u1%(u5%u8-u6%u7)-u4%(u2%u8-u4%u6)+u5%(u2%u7-u4%u5))+u3%(u1%(u4%u8-u6%u6)-u3%(u2%u8-u4%u6)+u5%(u2%u6-u4%u4))+u3%(u5%u8-u6%u7)-u4%(u4%u8-u6%u6)-u4%(u1%(u4%u7-u5%u6)-u3%(u2%u7-u4%u5)+u4%(u2%u6-u4%u4))+u5%(u4%u7-u5%u6))/det;
  vec b25= (u2%(u1%(u5%u7-u6%u6)-u4%(u2%u7-u3%u6)+u5%(u2%u6-u3%u5))-u3%(u1%(u4%u7-u5%u6)-u3%(u2%u7-u3%u6)+u5%(u2%u5-u3%u4))-u3%(u5%u7-u6%u6)+u4%(u4%u7-u5%u6)+u4%(u1%(u4%u6-u5%u5)-u3%(u2%u6-u3%u5)+u4%(u2%u5-u3%u4))-u5%(u4%u6-u5%u5))/det;
  
  vec b31= (u1%(u3%(u6%u8-u7%u7)-u5%(u4%u8-u5%u7)+u6%(u4%u7-u5%u6))-u2%(u2%(u6%u8-u7%u7)-u5%(u3%u8-u4%u7)+u6%(u3%u7-u4%u6))+u4%(u2%(u4%u8-u5%u7)-u3%(u3%u8-u4%u7)+(u3%u5-u4%u4)%u6)-u5%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u4%u6)+u5%(u3%u5-u4%u4)))/det;
  vec b32= (u1%(u2%(u6%u8-u7%u7)-u5%(u3%u8-u4%u7)+u6%(u3%u7-u4%u6))-u3%(u2%(u4%u8-u5%u7)-u3%(u3%u8-u4%u7)+(u3%u5-u4%u4)%u6)-u3%(u6%u8-u7%u7)+u5%(u4%u8-u5%u7)+u4%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u4%u6)+u5%(u3%u5-u4%u4))-u6%(u4%u7-u5%u6))/det;
  vec b33= (-u1%(u1%(u6%u8-u7%u7)-u4%(u3%u8-u4%u7)+u5%(u3%u7-u4%u6))+u3%(u1%(u4%u8-u5%u7)-u2%(u3%u8-u4%u7)+u5%(u3%u5-u4%u4))+u2%(u6%u8-u7%u7)-u4%(u4%u8-u5%u7)-u4%(u1%(u4%u7-u5%u6)-u2%(u3%u7-u4%u6)+u4%(u3%u5-u4%u4))+u5%(u4%u7-u5%u6))/det;
  vec b34= (u1%(u1%(u5%u8-u6%u7)-u4%(u2%u8-u4%u6)+u5%(u2%u7-u4%u5))-u3%(u1%(u3%u8-u5%u6)-u2%(u2%u8-u4%u6)+u5%(u2%u5-u3%u4))-u2%(u5%u8-u6%u7)+u4%(u3%u8-u5%u6)+u4%(u1%(u3%u7-u5%u5)-u2%(u2%u7-u4%u5)+u4%(u2%u5-u3%u4))-u5%(u3%u7-u5%u5))/det;
  vec b35= (-u1%(u1%(u5%u7-u6%u6)-u4%(u2%u7-u3%u6)+u5%(u2%u6-u3%u5))+u3%(u1%(u3%u7-u4%u6)-u2%(u2%u7-u3%u6)+(u2%u4-u3%u3)%u5)+u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)-u4%(u1%(u3%u6-u4%u5)-u2%(u2%u6-u3%u5)+u4%(u2%u4-u3%u3))+u5%(u3%u6-u4%u5))/det;
  
  vec b41= (-u1%(u3%(u5%u8-u6%u7)-u4%(u4%u8-u5%u7)+u6%(u4%u6-u5%u5))+u2%(u2%(u5%u8-u6%u7)-u4%(u3%u8-u4%u7)+u6%(u3%u6-u4%u5))-u3%(u2%(u4%u8-u5%u7)-u3%(u3%u8-u4%u7)+(u3%u5-u4%u4)%u6)+u5%(u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4)))/det;
  vec b42= (-u1%(u2%(u5%u8-u6%u7)-u4%(u3%u8-u4%u7)+u6%(u3%u6-u4%u5))+u2%(u2%(u4%u8-u5%u7)-u3%(u3%u8-u4%u7)+(u3%u5-u4%u4)%u6)+u3%(u5%u8-u6%u7)-u4%(u4%u8-u5%u7)-u4%(u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4))+u6%(u4%u6-u5%u5))/det;
  vec b43= (u1%(u1%(u5%u8-u6%u7)-u3%(u3%u8-u4%u7)+u5%(u3%u6-u4%u5))-u2%(u1%(u4%u8-u5%u7)-u2%(u3%u8-u4%u7)+u5%(u3%u5-u4%u4))-u2%(u5%u8-u6%u7)+u3%(u4%u8-u5%u7)+u4%(u1%(u4%u6-u5%u5)-u2%(u3%u6-u4%u5)+u3%(u3%u5-u4%u4))-u5%(u4%u6-u5%u5))/det;
  vec b44= (-u1%(u1%(u4%u8-u6%u6)-u3%(u2%u8-u4%u6)+u5%(u2%u6-u4%u4))+u2%(u1%(u3%u8-u5%u6)-u2%(u2%u8-u4%u6)+u5%(u2%u5-u3%u4))+u2%(u4%u8-u6%u6)-u3%(u3%u8-u5%u6)-u4%(u1%(u3%u6-u4%u5)-u2%(u2%u6-u4%u4)+u3%(u2%u5-u3%u4))+u5%(u3%u6-u4%u5))/det;
  vec b45= (u1%(u1%(u4%u7-u5%u6)-u3%(u2%u7-u3%u6)+u5%(u2%u5-u3%u4))-u2%(u1%(u3%u7-u4%u6)-u2%(u2%u7-u3%u6)+(u2%u4-u3%u3)%u5)-u2%(u4%u7-u5%u6)+u3%(u3%u7-u4%u6)+u4%(u1%(u3%u5-u4%u4)-u2%(u2%u5-u3%u4)+u3%(u2%u4-u3%u3))-u5%(u3%u5-u4%u4))/det;
  
  vec b51= (u1%(u3%(u5%u7-u6%u6)-u4%(u4%u7-u5%u6)+u5%(u4%u6-u5%u5))-u2%(u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)+u5%(u3%u6-u4%u5))+u3%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u4%u6)+u5%(u3%u5-u4%u4))-u4%(u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4)))/det;
  vec b52= (u1%(u2%(u5%u7-u6%u6)-u4%(u3%u7-u4%u6)+u5%(u3%u6-u4%u5))-u2%(u2%(u4%u7-u5%u6)-u3%(u3%u7-u4%u6)+u5%(u3%u5-u4%u4))-u3%(u5%u7-u6%u6)+u4%(u4%u7-u5%u6)+u3%(u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)+u4%(u3%u5-u4%u4))-u5%(u4%u6-u5%u5))/det;
  vec b53= (-u1%(u1%(u5%u7-u6%u6)-u3%(u3%u7-u4%u6)+u4%(u3%u6-u4%u5))+u2%(u1%(u4%u7-u5%u6)-u2%(u3%u7-u4%u6)+u4%(u3%u5-u4%u4))+u2%(u5%u7-u6%u6)-u3%(u4%u7-u5%u6)-u3%(u1%(u4%u6-u5%u5)-u2%(u3%u6-u4%u5)+u3%(u3%u5-u4%u4))+u4%(u4%u6-u5%u5))/det;
  vec b54= (u1%(u1%(u4%u7-u5%u6)-u3%(u2%u7-u4%u5)+u4%(u2%u6-u4%u4))-u2%(u1%(u3%u7-u5%u5)-u2%(u2%u7-u4%u5)+u4%(u2%u5-u3%u4))-u2%(u4%u7-u5%u6)+u3%(u3%u7-u5%u5)+u3%(u1%(u3%u6-u4%u5)-u2%(u2%u6-u4%u4)+u3%(u2%u5-u3%u4))-u4%(u3%u6-u4%u5))/det;
  vec b55= (-u1%(u1%(u4%u6-u5%u5)-u3%(u2%u6-u3%u5)+u4%(u2%u5-u3%u4))+u2%(u1%(u3%u6-u4%u5)-u2%(u2%u6-u3%u5)+u4%(u2%u4-u3%u3))+u2%(u4%u6-u5%u5)-u3%(u3%u6-u4%u5)-u3%(u1%(u3%u5-u4%u4)-u2%(u2%u5-u3%u4)+u3%(u2%u4-u3%u3))+u4%(u3%u5-u4%u4))/det;
  '
  
  
  switch(Dindex,
{
  Dpart='
  vec p=(1.0/3.0) *(3*(x0.col(3)/6.0)%x0.col(1) - pow(x0.col(2)/2.0,2))/pow(x0.col(3)/6.0,2);
  vec q=(1.0/27.0)*(27*pow(x0.col(3)/6.0,2)%(x0.col(0)-Xt) - 9*(x0.col(3)/6.0)%(x0.col(2)/2.0)%x0.col(1) + 2*pow(x0.col(2)/2.0,3))/pow(x0.col(3)/6.0,3);
  vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  vec th=-(x0.col(2)/2.0)/(3*(x0.col(3)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));
  
  vec K =x0.col(0)%th+(x0.col(1)%th%th)/2.0+(x0.col(2)%th%th%th)/6.0 +(x0.col(3)%th%th%th%th)/24.0;
  vec K1=x0.col(0)   +(x0.col(1)%th)       +(x0.col(2)%th%th)/2.0    +(x0.col(3)%th%th%th)/6.0;
  vec K2=x0.col(1)   +(x0.col(2)%th)       +(x0.col(3)%th%th)/2.0;
  vec val=-0.5*log(2*3.141592653589793*K2)+(K-th%K1);
  return(val);
}
  '
},
{
  Dpart=
    '
  vec betas1 =+b12+2*b13%u1+3*b14%u2+4*b15%u3;
  vec betas2 =+b22+2*b23%u1+3*b24%u2+4*b25%u3;
  vec betas3 =+b32+2*b33%u1+3*b34%u2+4*b35%u3;
  vec betas4 =+b42+2*b43%u1+3*b44%u2+4*b45%u3;
  vec betas5 =+b52+2*b53%u1+3*b54%u2+4*b55%u3;
  
  vec lo =-(0.5/(u1-lower))%(sqrt(exp(2*alpha)+4*pow((u1-lower),2))-exp(alpha));
  vec up =-(0.5/(u1-upper))%(sqrt(exp(2*alpha)+4*pow((u1-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  vec K=exp(-betas1%pow(tau,1)-0.5*betas2%pow(tau,2)-0.333333333333333333*betas3%pow(tau,3)-0.25*betas4%pow(tau,4)-0.2*betas5%pow(tau,5))%rho%DT;
  for (int i = 1; i <= P; i++)
{
  lo=lo+DT;
  tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  K=K+exp(-betas1%pow(tau,1)-0.5*betas2%pow(tau,2)-0.333333333333333333*betas3%pow(tau,3)-0.25*betas4%pow(tau,4)-0.2*betas5%pow(tau,5))%rho%DT;
}
  
  return((-betas1%pow(Xt,1)-0.5*betas2%pow(Xt,2)-0.333333333333333333*betas3%pow(Xt,3)-0.25*betas4%pow(Xt,4)-0.2*betas5%pow(Xt,5))-log(K));
}
  '
},
{
  Dpart='
  vec betas1 =+b11+2*b12%u1+3*b13%u2+4*b14%u3+4*b15%u4;
  vec betas2 =+b21+2*b22%u1+3*b23%u2+4*b24%u3+4*b25%u4;
  vec betas3 =+b31+2*b32%u1+3*b33%u2+4*b34%u3+4*b35%u4;
  vec betas4 =+b41+2*b42%u1+3*b43%u2+4*b44%u3+4*b45%u4;
  vec betas5 =+b51+2*b52%u1+3*b53%u2+4*b54%u3+4*b55%u4;
  
  vec lo =-(0.5/(u1-lower))%(sqrt(exp(2*alpha)+4*pow((u1-lower),2))-exp(alpha));
  vec up =-(0.5/(u1-upper))%(sqrt(exp(2*alpha)+4*pow((u1-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  
  vec K=exp(-betas1%log(tau)+(-betas2%tau-0.5*betas3%pow(tau,2)-0.33333333333333333333333*betas4%pow(tau,3)-0.25*betas5%pow(tau,4)))%rho%DT;
  for (int i = 1; i <= P; i++)
{
  lo=lo+DT;
  tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  K=K+exp(-betas1%log(tau)+(-betas2%tau-0.5*betas3%pow(tau,2)-0.33333333333333333333333*betas4%pow(tau,3)-0.25*betas5%pow(tau,4)))%rho%DT;
}
  return((-betas1%log(Xt)+(-betas2%Xt-0.5*betas3%pow(Xt,2)-0.33333333333333333333333*betas4%pow(Xt,3)-0.25*betas5%pow(Xt,4)))-log(K));
}
  '
  },
{
  Dpart='
  vec betas1 =+2*b11%u1+3*b12%u2+4*b13%u3+5*b14%u4+6*b15%u5;
  vec betas2 =+2*b21%u1+3*b22%u2+4*b23%u3+5*b24%u4+6*b25%u5;
  vec betas3 =+2*b31%u1+3*b32%u2+4*b33%u3+5*b34%u4+6*b35%u5;
  vec betas4 =+2*b41%u1+3*b42%u2+4*b43%u3+5*b44%u4+6*b45%u5;
  vec betas5 =+2*b51%u1+3*b52%u2+4*b53%u3+5*b54%u4+6*b55%u5;
  
  vec lo =-(0.5/(u1-lower))%(sqrt(exp(2*alpha)+4*pow((u1-lower),2))-exp(alpha));
  vec up =-(0.5/(u1-upper))%(sqrt(exp(2*alpha)+4*pow((u1-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  
  vec K=exp(-betas2%log(tau)+(betas1/tau-betas3%tau-0.5*betas4%pow(tau,2)-0.3333333333333333*betas5%pow(tau,3)))%rho%DT;
  for (int i = 1; i <= P; i++)
{
  lo=lo+DT;
  tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  K=K+exp(-betas2%log(tau)+(betas1/tau-betas3%tau-0.5*betas4%pow(tau,2)-0.3333333333333333*betas5%pow(tau,3)))%rho%DT;
}
  return((-betas2%log(Xt)+(betas1/Xt-betas3%Xt-0.5*betas4%pow(Xt,2)-0.3333333333333333*betas5%pow(Xt,3)))-log(K));
}
  '
  },
{
  Dpart='
  
  vec betas1 =+b11%(1-2*u1)+b12%(2*u1-3*u2)+b13%(3*u2-4*u3)+b14%(4*u3-5*u4)+b15%(5*u4-6*u5);
  vec betas2 =+b21%(1-2*u1)+b22%(2*u1-3*u2)+b23%(3*u2-4*u3)+b24%(4*u3-5*u4)+b25%(5*u4-6*u5);
  vec betas3 =+b31%(1-2*u1)+b32%(2*u1-3*u2)+b33%(3*u2-4*u3)+b34%(4*u3-5*u4)+b35%(5*u4-6*u5);
  vec betas4 =+b41%(1-2*u1)+b42%(2*u1-3*u2)+b43%(3*u2-4*u3)+b44%(4*u3-5*u4)+b45%(5*u4-6*u5);
  vec betas5 =+b51%(1-2*u1)+b52%(2*u1-3*u2)+b53%(3*u2-4*u3)+b54%(4*u3-5*u4)+b55%(5*u4-6*u5);
  
  vec lo =-(0.5/(u1-lower))%(sqrt(exp(2*alpha)+4*pow((u1-lower),2))-exp(alpha));
  vec up =-(0.5/(u1-upper))%(sqrt(exp(2*alpha)+4*pow((u1-upper),2))-exp(alpha));
  vec DT = (up-lo)/(1.00*P);
  vec tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  vec rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  
  vec K=exp(-betas2%log(tau)+(betas1/tau-betas3%tau-0.5*betas4%pow(tau,2)-0.3333333333333333*betas5%pow(tau,3)))%rho%DT;
  for (int i = 1; i <= P; i++)
{
  lo=lo+DT;
  tau   = exp(alpha)*lo/(1-pow(lo,2))+u1;
  rho   = exp(alpha)*(1+pow(lo,2))/pow((1-pow(lo,2)),2);
  K=K+exp(-betas2%log(tau)+(betas1/tau-betas3%tau-0.5*betas4%pow(tau,2)-0.3333333333333333*betas5%pow(tau,3)))%rho%DT;
}
  return((-betas2%log(Xt)+(betas1/Xt-betas3%Xt-0.5*betas4%pow(Xt,2)-0.3333333333333333*betas5%pow(Xt,3)))-log(K));
}
  '
  })
if(Dindex!=1)
{
  Dpart=paste(Inv,Dpart)
}
}


#==============================================================================
#                                Syntax Check
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
    stop('^-Operator not defined in C++: Use pow(x,p) instead (e.g. pow(x,2) i.s.o. x^2).',call. = F)
    k=k+1
  }
  return(k)
}
for(i in which(func.list==1))
{
  strcheck(body(namess[i])[2])
}
#==============================================================================
#                   Generate TYPE of Solution
#==============================================================================
if(state1)
{
  # DATA RESOLUTION -------------------------------------------------------------
  if(!homo.res)
  {
    delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh)
    if(!is.null(exclude))
    {
      diffs=diff(T.seq)/mesh
      diffs[excl]=1/20/mesh
      delt=cbind(diffs,diffs)
    }
  }
  # REFERENCE MATRIX ------------------------------------------------------------
  MAT=rbind(
    c('(1+0*a.col(0))','a.col(0)' ,''),
    c('','2*a.col(1)','(1+0*a.col(0))'),
    c('','3*a.col(2)',''),
    c('','4*a.col(3)',''))
  
  # HOMOGENEITY -----------------------------------------------------------------
  namess2=c('G0','G1','Q0')
  func.list2=rep(0,3)
  obs=objects(pos=1)
  for(i in 1:3)
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
  
  #func.list.timehomo[c(3)]=1    # Always set these to 1
  
  #BUILD ODE --------------------------------------------------------------------
  dims=rep('(',2)
  for(i in 1:2)
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
  
  if(any(dims=='()')){dims[which(dims=='()')]='0*a.col(0)'}
  
  
  for(i in 1:2)
  {
    dims[i]=paste0(paste0(paste0('   atemp.col(',i-1,')='),dims[i]),';')
  }
  # WRIGHT AND SOURCE -----------------------------------------------------------
  
  txt.full=paste(fpart,'\n',dims[1],'\n',dims[2],ODEpart,Dpart)
  type.sol ="                 Generalized Ornstein-Uhlenbeck "
  #library(Rcpp)
  #library(RcppArmadillo)
  if(wrt){write(txt.full,file='GQD.mcmc.cpp')}
  stre="Compiling C++ code. Please wait."
  cat(stre, " \r")
  flush.console()
  sourceCpp(code=txt.full)
  cat('                                     ','\r')
}
if(state2)
{
  # DATA RESOLUTION -------------------------------------------------------------
  if((!homo.res)&(TR.order==4))
  {
    delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh)
    if(!is.null(exclude))
    {
      diffs=diff(T.seq)/mesh
      diffs[excl]=1/20/mesh
      delt=cbind(diffs,diffs,diffs,diffs)
    }
  }
  if((!homo.res)&(TR.order==6))
  {
    delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh)
    if(!is.null(exclude))
    {
      diffs=diff(T.seq)/mesh
      diffs[excl]=1/20/mesh
      delt=cbind(diffs,diffs,diffs,diffs,diffs,diffs)
    }
  }
  if((!homo.res)&(TR.order==8))
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
  if(TR.order==4)
  {
    MAT=rbind(
      c('(1+0*a.col(0))','a.col(0)','(a.col(1)+a.col(0)%a.col(0))','','',''),
      c('','(2*a.col(1))','(2*a.col(2)+4*a.col(0)%a.col(1))','(1+0*a.col(0))','a.col(0)','(a.col(1)+a.col(0)%a.col(0))'),
      c('','(3*a.col(2))','(3*a.col(3)+6*a.col(0)%a.col(2)+6*a.col(1)%a.col(1))','','(3*a.col(1))','(3*a.col(2)+6*a.col(0)%a.col(1))'),
      c('','(4*a.col(3))','(8*a.col(0)%a.col(3)+24*a.col(1)%a.col(2))','','(6*a.col(2))','(6*a.col(3)+12*a.col(0)%a.col(2)+12*a.col(1)%a.col(1))'))
  }
  if(TR.order==6)
  {
    MAT = rbind(
      c('(1+0*a.col(0))','(a.col(0))'  ,'(1*a.col(1)+1*a.col(0)%a.col(0))'                                            ,'','',''),
      c('','(2*a.col(1))','(2*a.col(2)+4*a.col(0)%a.col(1))'                                            ,'(1+0*a.col(0))','(   a.col(0))','(a.col(1)+a.col(0)%a.col(0))'),
      c('','(3*a.col(2))','(3*a.col(3)+6*a.col(0)%a.col(2)+6*a.col(1)%a.col(1))'                        ,'','( 3*a.col(1))','(3*a.col(2)+6*a.col(0)%a.col(1))'),
      c('','(4*a.col(3))','(4*a.col(4)+8*a.col(0)%a.col(3)+24*a.col(1)%a.col(2))'                       ,'','( 6*a.col(2))','(6*a.col(3)+12*a.col(0)%a.col(2)+12*a.col(1)%a.col(1))'),
      c('','(5*a.col(4))','(5*a.col(5)+10*a.col(0)%a.col(4)+40*a.col(1)%a.col(3)+30*a.col(2)%a.col(2))' ,'','(10*a.col(3))','(10*a.col(4)+20*a.col(0)%a.col(3)+60*a.col(1)%a.col(2))'),
      c('','(6*a.col(5))','(          +12*a.col(0)%a.col(5)+60*a.col(1)%a.col(4)+120*a.col(2)%a.col(3))','','(15*a.col(4))','(15*a.col(5)+30*a.col(0)%a.col(4)+120*a.col(1)%a.col(3)+90*a.col(2)%a.col(2))'))
  }
  if(TR.order==8)
  {
    MAT = rbind(
      c('(1+0*a.col(0))','  (a.col(0))','(1*a.col(1) +1*a.col(0)%a.col(0))'                                                                          ,'','',''),
      c('','(2*a.col(1))','(2*a.col(2) +4*a.col(0)%a.col(1))'                                                                          ,'(1+0*a.col(0))','(   a.col(0))','(a.col(1)+a.col(0)%a.col(0))'),
      c('','(3*a.col(2))','(3*a.col(3) +6*a.col(0)%a.col(2)  +6*a.col(1)%a.col(1))'                                                  ,'','( 3*a.col(1))','( 3*a.col(2) +6*a.col(0)%a.col(1))'),
      c('','(4*a.col(3))','(4*a.col(4) +8*a.col(0)%a.col(3) +24*a.col(1)%a.col(2))'                                                  ,'','( 6*a.col(2))','( 6*a.col(3)+12*a.col(0)%a.col(2) +12*a.col(1)%a.col(1))'),
      c('','(5*a.col(4))','(5*a.col(5)+10*a.col(0)%a.col(4) +40*a.col(1)%a.col(3)   +30*a.col(2)%a.col(2))'                        ,'','(10*a.col(3))','(10*a.col(4)+20*a.col(0)%a.col(3) +60*a.col(1)%a.col(2))'),
      c('','(6*a.col(5))','(6*a.col(6)+12*a.col(0)%a.col(5) +60*a.col(1)%a.col(4)  +120*a.col(2)%a.col(3))'                        ,'','(15*a.col(4))','(15*a.col(5)+30*a.col(0)%a.col(4)+120*a.col(1)%a.col(3) +90*a.col(2)%a.col(2))'),
      c('','(7*a.col(6))','(7*a.col(7)+14*a.col(0)%a.col(6) +84*a.col(1)%a.col(5)  +210*a.col(2)%a.col(4)+140*a.col(3)%a.col(3))','','(21*a.col(5))','(21*a.col(6)+42*a.col(0)%a.col(5)+210*a.col(1)%a.col(4)+420*a.col(2)%a.col(3))'),
      c('','(8*a.col(7))','(           +16*a.col(0)%a.col(7)+112*a.col(1)%a.col(6)  +336*a.col(2)%a.col(5)+560*a.col(4)%a.col(3))','','(28*a.col(6))','(28*a.col(7)+56*a.col(0)%a.col(6)+336*a.col(1)%a.col(5)+840*a.col(2)%a.col(4)+560*a.col(3)%a.col(3))'))
  }
  
  # HOMOGENEITY -----------------------------------------------------------------
  namess2=c('G0','G1','G2','Q0','Q1','Q2')
  func.list2=rep(0,6)
  obs=objects(pos=1)
  for(i in 1:6)
  {
    if(sum(obs==namess[i])){func.list2[i]=1}
  }
  
  func.list.timehomo=func.list2*0
  
  for(i in which(func.list2==1))
  {
    # which expressions vary over time
    result=eval(body(namess2[i]))
    func.list.timehomo[i]=2-(sum(diff(result)==0)==(length(result)-1))
  }
  
  if(any(func.list.timehomo==2)){homo=F}
  
  #func.list.timehomo[c(4)]=1    # Always set these to 1
  
  #BUILD ODE --------------------------------------------------------------------
  dims=rep('(',TR.order)
  for(i in 1:TR.order)
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
  
  if(any(dims=='()')){dims[which(dims=='()')]='0*a.col(0)'}
  
  
  for(i in 1:TR.order)
  {
    dims[i]=paste0(paste0(paste0('   atemp.col(',i-1,')='),dims[i]),';')
  }
  # WRIGHT AND SOURCE -----------------------------------------------------------
  
  
  if(TR.order==4)
  {
    txt.full=paste(fpart,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],ODEpart,Dpart)
  }
  if(TR.order==6)
  {
    txt.full=paste(fpart,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],'\n',dims[6],ODEpart,Dpart)
  }
  if(TR.order==8)
  {
    txt.full=paste(fpart,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],'\n',dims[6],'\n',dims[7],'\n',dims[8],ODEpart,Dpart)
  }
  type.sol ="                 Generalized Quadratic Diffusion (GQD) "
  
  #library(Rcpp)
  #library(RcppArmadillo)
  if(wrt){write(txt.full,file='GQD.mcmc.cpp')}
  stre="Compiling C++ code. Please wait."
  cat(stre, " \r")
  flush.console()
  sourceCpp(code=txt.full)
  cat('                                     ','\r')
}


   #==============================================================================
   #                           Interface Module
   #==============================================================================

    namess4=namess
    trim <- function (x) gsub("([[:space:]])", "", x)

    for(i in 1:6)
    {
        if(sum(obs==namess4[i]))
        {
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
    buffer7=c('_______________________ Model/Chain Info _______________________')

    Info=c(buffer0,type.sol,buffer0,buffer4,namess4[1:3],buffer5,namess4[4:6])
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
like=function(pars)
{
  sum(solver(X[-nnn],X[-1],c(0,pars),mesh,delt,nnn-1,T.seq[-nnn],P,alpha,lower,upper,TR.order)[-excl])
}

    tme=Sys.time()
    result=optim(theta,like,control = control,method=method, hessian = T)
    tme=tme.eval(tme)

    actual.p=length(theta)
    model.inf=list(elapsed.time=tme,time.homogeneous=c('Yes','No')[2-homo],p=actual.p,N=length(X)-length(excl)+1,Tag=Tag)
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
    class(ret) = 'GQD.mle'

    return(ret)
  }
