globalVariables(c('priors','Lam0','Lam1','Jbeta','Jmu','Jsig','Jalpha','Jlam','Ja','Jb'))
JGQD.mcmc <-
function(X,time,mesh=10,theta,sds,updates=1000,burns=min(round(updates/2),25000),Jtype='Add',Jdist='Normal',Dtype='Saddlepoint',RK.order=4,exclude=NULL,plot.chain=TRUE,wrt=FALSE,Tag=NA,factorize=TRUE,print.output=TRUE)
{
  P=100;alpha=0.1;lower=0;upper=50
  Trunc=c(4,4)
  solver   =function(Xs, Xt, theta, N , delt , N2, tt  , P , alpha, lower , upper, tro){}
  rm(list =c('solver'))

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
   Dtypes =c('Saddlepoint')
   JDtypes=c('Normal','Exponential','Gamma','Laplace')
   Dindex = which(Dtypes==Dtype)
   JDindex = which(JDtypes==Jdist)
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
    ,'18. Density: Dtype has to be of type Saddlepoint.\n'
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
    ,'39. Input: {Jdist} has to be of type Normal, Exponential, Gamma or Laplace.\n'
    ,'40. Input: {Jtype} has to be of type Add or Mult.\n'
    ,'41. Input: {factorize} has to be TRUE or FALSE.\n'
    ,'42. Input: Current version supports {Trunc[1]} = 4 or 8.\n'
    ,'43. Input: Current version supports {Trunc[2]} = 4.\n'
  )

   warntrue = rep(F,50)
      check.thetas = function(theta,tt)
   {
     t=tt
     theta = theta+runif(length(theta),0.001,0.002)*sign(theta)
     namess=c('G0','G1','G2','Q0','Q1','Q2','Lam0','Lam1','Jmu','Jsig','Jlam','Jalpha','Jbeta','Ja','Jb')
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
     namess=c('G0','G1','G2','Q0','Q1','Q2','Lam0','Lam1','Jmu','Jsig','Jlam','Jalpha','Jbeta','Ja','Jb')
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
   if(missing(X))                                                {warntrue[1]=TRUE}
   if(missing(time))                                             {warntrue[2]=TRUE}
   if(missing(theta))                                            {warntrue[3]=TRUE}
   if(missing(sds))                                              {warntrue[4]=TRUE}
   if(!is.vector(X))                                             {warntrue[5]=TRUE}
   if(!is.vector(time))                                          {warntrue[6]=TRUE}
   # Check model parameters:
   if(check.thetas2(theta)!=0)                                   {warntrue[7]=TRUE}
   if(!warntrue[7]){if(any(check.thetas(theta,T.seq)==0))        {warntrue[8]=TRUE}}

   # Check input length:
   if(length(X)<10)                                              {warntrue[9]=TRUE}
   if(length(time)<10)                                          {warntrue[10]=TRUE}
   if(length(lower)>1)                                          {warntrue[11]=TRUE}
   if(length(upper)>1)                                          {warntrue[12]=TRUE}
   if(length(P)!=1)                                             {warntrue[13]=TRUE}
   if(length(mesh)!=1)                                        {warntrue[14]=TRUE}
   if(length(alpha)!=1)                                         {warntrue[15]=TRUE}
   if(length(Trunc)!=2)                                         {warntrue[16]=TRUE}
   if(length(RK.order)!=1)                                      {warntrue[17]=TRUE}
   if(length(updates)!=1)                                     {warntrue[33]=TRUE}
   if(length(burns)!=1)                                         {warntrue[34]=TRUE}
   if(length(Dtype)!=1)                                         {warntrue[37]=TRUE}


   # Check density approx parameters:
   if(sum(Dindex)==0)                                          {warntrue[18] =TRUE}
   if(!warntrue[18])
   {
    if((Dindex==3)|(Dindex==4)){if(lower[1]<=0)                {warntrue[19] =TRUE}}
    if(Dindex==5){if(any(X<=0)|any(X>=1))                      {warntrue[20] =TRUE}}
   }
   if(!any(warntrue[c(11,12)])){if(upper<=lower)                {warntrue[21] =TRUE}}
   if(!warntrue[13]){if(P<10)                                   {warntrue[22] =TRUE}}
   if(!warntrue[16]){if(Trunc[2]>Trunc[1])                      {warntrue[23] =TRUE}}
   #if(sum(c(4,8)==Trunc[1])==0)                                 {warntrue[42] =TRUE}}
   #if(sum(c(4)==Trunc[2])  ==0)                                 {warntrue[43] =TRUE}}
   #  Miscelaneous checks:
   excl=0
   if(is.null(exclude)){excl=length(T.seq)-1+200}
   if(!is.null(exclude)){excl=exclude}
   test.this =max(diff(T.seq)[-excl])/mesh
   if(test.this>0.1)                                            {warntrue[24]=TRUE}
   if(test.this>=1)                                             {warntrue[25]=TRUE}
   if(!warntrue[17]){if(!((RK.order==4)|(RK.order==10)))        {warntrue[26]=TRUE}}
   if(!warntrue[14]){if(mesh<5)                               {warntrue[27]=TRUE}}
   if(length(X)!=length(time))                                  {warntrue[28]=TRUE}
   if(!any(warntrue[c(33,34)])){if(burns>updates)             {warntrue[29]=TRUE}}
   if(!warntrue[33]){if(updates<2)                            {warntrue[30]=TRUE}}
   if(length(theta)!=length(sds))                               {warntrue[31]=TRUE}
   if(any(is.na(X))||any(is.na(time)))                          {warntrue[36]=TRUE}
   if(sum(JDindex)==0)                                          {warntrue[39] =TRUE}
   if((Jtype!='Add')&&(Jtype!='Mult'))                          {warntrue[40] =TRUE}
   if((factorize!=TRUE)&&(factorize!=FALSE))                    {warntrue[41] =TRUE}
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


   nnn=length(X)

   #==============================================================================
   #                           Data resolution Module
   #==============================================================================
   #t=seq(0,100,by=1/100)
   homo=TRUE
   homo.res=TRUE
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
     prior.list=paste0('d(theta)',':',' No priors given.')
     pp=function(theta){1}
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
   if(state1){DTR.order=2;TR.order=4;sol.state='Normally distributed diffusion.';}
   if((state1&(Dtype!='Saddle'))){TR.order=4;DTR.order=2;sol.state='2nd Ord. Truncation + Std Normal Dist.';}
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
    // [[Rcpp::export]
    mat f(mat a,vec theta,vec t,int N2)
    {

       mat atemp(N2,10);'
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

       mat atemp(N2,18);'
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
  mat fx0(N2,2*tro+2);
  mat fx1(N2,2*tro+2);
  mat fx2(N2,2*tro+2);
  mat fx3(N2,2*tro+2);
  mat fx4(N2,2*tro+2);
  mat fx5(N2,2*tro+2);
  mat fx6(N2,2*tro+2);
  mat fx7(N2,2*tro+2);
  mat fx8(N2,2*tro+2);
  mat fx9(N2,2*tro+2);
  mat fx10(N2,2*tro+2);
  mat fx11(N2,2*tro+2);
  mat fx12(N2,2*tro+2);
  mat fx13(N2,2*tro+2);
  mat fx14(N2,2*tro+2);
  mat fx15(N2,2*tro+2);
  mat fx16(N2,2*tro+2);
  mat x0(N2,2*tro+2);
  mat y0(N2,2*tro+2);
  mat x1(N2,2*tro+2);
  mat x2(N2,2*tro+2);
  mat x3(N2,2*tro+2);
  mat x4(N2,2*tro+2);
  mat x5(N2,2*tro+2);
  mat x6(N2,2*tro+2);
  mat x7(N2,2*tro+2);
  mat x8(N2,2*tro+2);
  mat x9(N2,2*tro+2);
  mat x10(N2,2*tro+2);
  mat x11(N2,2*tro+2);
  mat x12(N2,2*tro+2);
  mat x13(N2,2*tro+2);
  mat x14(N2,2*tro+2);
  mat x15(N2,2*tro+2);
  mat x16(N2,2*tro+2);
  x0.fill(0);
  for (int i = 1; i < tro+1; i++)
  {
    x0.col(i-1)=pow(Xs,i);
    x0.col(i-1+tro)=pow(Xs,i);
  }
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
  mat fx0(N2,2*tro+2);
  mat fx1(N2,2*tro+2);
  mat fx2(N2,2*tro+2);
  mat fx3(N2,2*tro+2);
  mat fx4(N2,2*tro+2);
  mat fx5(N2,2*tro+2);
  mat fx6(N2,2*tro+2);
  mat fx7(N2,2*tro+2);
  mat fx8(N2,2*tro+2);
  mat fx9(N2,2*tro+2);
  mat fx10(N2,2*tro+2);
  mat fx11(N2,2*tro+2);
  mat fx12(N2,2*tro+2);
  mat fx13(N2,2*tro+2);
  mat fx14(N2,2*tro+2);
  mat fx15(N2,2*tro+2);
  mat fx16(N2,2*tro+2);
  mat x0(N2,2*tro+2);
  mat y0(N2,2*tro+2);
  mat x1(N2,2*tro+2);
  mat x2(N2,2*tro+2);
  mat x3(N2,2*tro+2);
  mat x4(N2,2*tro+2);
  mat x5(N2,2*tro+2);
  mat x6(N2,2*tro+2);
  mat x7(N2,2*tro+2);
  mat x8(N2,2*tro+2);
  mat x9(N2,2*tro+2);
  mat x10(N2,2*tro+2);
  mat x11(N2,2*tro+2);
  mat x12(N2,2*tro+2);
  mat x13(N2,2*tro+2);
  mat x14(N2,2*tro+2);
  mat x15(N2,2*tro+2);
  mat x16(N2,2*tro+2);
  x0.fill(0);
  for (int i = 1; i < tro+1; i++)
  {
    x0.col(i-1)=pow(Xs,i);
    x0.col(i-1+tro)=pow(Xs,i);
  }
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
List  solver(vec Xs,vec Xt,vec theta,int N,double delt,int N2,vec tt,int P,double alpha,double lower,double upper,int tro,double eps=0.01)
{
    int steps=0;
    mat x0(N2,2*tro+2);
  mat y0(N2,2*tro+2);
  mat xa(N2,2*tro+2);
  mat xb(N2,2*tro+2);
  mat fx1(N2,2*tro+2);
  mat fx2(N2,2*tro+2);
  mat fx3(N2,2*tro+2);
  mat fx4(N2,2*tro+2);
  mat fx5(N2,2*tro+2);
  mat fx6(N2,2*tro+2);

  mat etemp(N2,1);
  x0.fill(0);
  for (int i = 1; i < tro+1; i++)
  {
  x0.col(i-1)=pow(Xs,i);
  x0.col(i-1+tro)=pow(Xs,i);
  }
  vec dd = tt;
  double d=tt(0);
  double t=tt(1);

  for (int i = 1; i < N+1; i++)
  {
  delt = std::min(delt,t-d);
  fx1 = delt*f(x0,theta,dd,N2);
  fx2 = delt*f(x0+0.25*fx1,theta,dd+0.25*delt,N2);
  fx3 = delt*f(x0+0.09375*fx1+0.28125*fx2,theta,dd+0.375*delt,N2);
  fx4 = delt*f(x0+0.879381*fx1-3.277196*fx2+ 3.320892*fx3,theta,dd+0.9230769*delt,N2);
  fx5 = delt*f(x0+2.032407*fx1-8*fx2+7.173489*fx3-0.2058967*fx4,theta,dd+delt,N2);
  fx6 = delt*f(x0-0.2962963*fx1+2*fx2-1.381676*fx3+0.4529727*fx4-0.275*fx5,theta,dd+0.5*delt,N2);
  xa=  x0+0.1157407*fx1+0.5489279*fx3+0.5353314*fx4-0.2*fx5;
  xb =  x0+0.1185185*fx1+0.5189864*fx3+0.5061315*fx4-0.18*fx5+0.03636364*fx6;

  x0=xa;
  d=d+delt;
  dd=dd+delt;
  }
    '
 }
 if(!homo.res)
 {
 ODEpart= '
    return atemp;
}

// [[Rcpp::export]]
List  solver(vec Xs,vec Xt,vec theta,int N,mat delt,int N2,vec tt,int P,double alpha,double lower,double upper,int tro,double eps=0.01)
{
    int steps=0;
    mat x0(N2,2*tro+2);
  mat y0(N2,2*tro+2);
  mat xa(N2,2*tro+2);
  mat xb(N2,2*tro+2);
  mat fx1(N2,2*tro+2);
  mat fx2(N2,2*tro+2);
  mat fx3(N2,2*tro+2);
  mat fx4(N2,2*tro+2);
  mat fx5(N2,2*tro+2);
  mat fx6(N2,2*tro+2);

  mat etemp(N2,1);
  x0.fill(0);
  for (int i = 1; i < tro+1; i++)
  {
  x0.col(i-1)=pow(Xs,i);
  x0.col(i-1+tro)=pow(Xs,i);
  }
  vec dd = tt;

  for (int i = 1; i < N+1; i++)
  {

  fx1 = delt%f(x0,theta,dd,N2);
  fx2 = delt%f(x0+0.25*fx1,theta,dd+0.25*delt.col(0),N2);
  fx3 = delt%f(x0+0.09375*fx1+0.28125*fx2,theta,dd+0.375*delt.col(0),N2);
  fx4 = delt%f(x0+0.879381*fx1-3.277196*fx2+ 3.320892*fx3,theta,dd+0.9230769*delt.col(0),N2);
  fx5 = delt%f(x0+2.032407*fx1-8*fx2+7.173489*fx3-0.2058967*fx4,theta,dd+delt.col(0),N2);
  fx6 = delt%f(x0-0.2962963*fx1+2*fx2-1.381676*fx3+0.4529727*fx4-0.275*fx5,theta,dd+0.5*delt.col(0),N2);
  xa=  x0+0.1157407*fx1+0.5489279*fx3+0.5353314*fx4-0.2*fx5;
  xb =  x0+0.1185185*fx1+0.5189864*fx3+0.5061315*fx4-0.18*fx5+0.03636364*fx6;
  x0=xa;
  dd=dd+delt.col(0);
  }
    '
 }
}

if(DTR.order==2)
{
    Dpart='
    y0.col(0)= x0.col(0+4);
    y0.col(1)= x0.col(1+4)-1*y0.col(0)%x0.col(0+4);
    y0.col(2)= x0.col(2+4)-1*y0.col(0)%x0.col(1+4)-2*y0.col(1)%x0.col(0+4);
    y0.col(3)= x0.col(3+4)-1*y0.col(0)%x0.col(2+4)-3*y0.col(1)%x0.col(1+4)-3*y0.col(2)%x0.col(0+4);

    y0.col(4)= ((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
    y0.col(5)= ((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
    y0.col(6)= ((x0.col(2)-x0.col(6)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-2*y0.col(5)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
    y0.col(7)= ((x0.col(3)-x0.col(7)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(2)-x0.col(6)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-3*y0.col(5)%((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-3*y0.col(6)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));

    vec val=exp(-0.5*log(2*3.141592653589793*y0.col(1))-0.5*((Xt-y0.col(0))%(Xt-y0.col(0))/y0.col(1))-x0.col(8)-x0.col(9));
    vec p=(1.0/3.0) *(3*(y0.col(7)/6.0)%y0.col(5) - pow(y0.col(6)/2.0,2))/pow(y0.col(7)/6.0,2);
    vec q=(1.0/27.0)*(27*pow(y0.col(7)/6.0,2)%(y0.col(4)-Xt) - 9*(y0.col(7)/6.0)%(y0.col(6)/2.0)%y0.col(5) + 2*pow(y0.col(6)/2.0,3))/pow(y0.col(7)/6.0,3);
    vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
    vec th=-(y0.col(6)/2.0)/(3*(y0.col(7)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

    vec  K =y0.col(4)%th+(y0.col(5)%th%th)/2.0+(y0.col(6)%th%th%th)/6.0 +(y0.col(7)%th%th%th%th)/24.0;
    vec  K1=y0.col(4)   +(y0.col(5)%th)       +(y0.col(6)%th%th)/2.0    +(y0.col(7)%th%th%th)/6.0;
    vec  K2=y0.col(5)   +(y0.col(6)%th)       +(y0.col(7)%th%th)/2.0;
    val= log(val+exp(-0.5*log(2*3.141592653589793*K2)+(K-th%K1))%(1-exp(-x0.col(8)-x0.col(9))));

    List ret;
    ret["like"]  = val;
    ret["exc"]  = (1-exp(-x0.col(8)-x0.col(9)));
    ret["steps"]  = steps;
    return(ret);
}'
}
if(DTR.order==4)
{
  if(factorize)
  {
  switch(Dindex,
{
  Dpart='
  y0.col(0)= x0.col(0+4);
  y0.col(1)= x0.col(1+4)-1*y0.col(0)%x0.col(0+4);
  y0.col(2)= x0.col(2+4)-1*y0.col(0)%x0.col(1+4)-2*y0.col(1)%x0.col(0+4);
  y0.col(3)= x0.col(3+4)-1*y0.col(0)%x0.col(2+4)-3*y0.col(1)%x0.col(1+4)-3*y0.col(2)%x0.col(0+4);

  y0.col(4)= ((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
  y0.col(5)= ((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
  y0.col(6)= ((x0.col(2)-x0.col(6)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-2*y0.col(5)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
  y0.col(7)= ((x0.col(3)-x0.col(7)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(2)-x0.col(6)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-3*y0.col(5)%((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-3*y0.col(6)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));

  vec p=(1.0/3.0) *(3*(y0.col(3)/6.0)%y0.col(1) - pow(y0.col(2)/2.0,2))/pow(y0.col(3)/6.0,2);
  vec q=(1.0/27.0)*(27*pow(y0.col(3)/6.0,2)%(y0.col(0)-Xt) - 9*(y0.col(3)/6.0)%(y0.col(2)/2.0)%y0.col(1) + 2*pow(y0.col(2)/2.0,3))/pow(y0.col(3)/6.0,3);
  vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  vec th=-(y0.col(2)/2.0)/(3*(y0.col(3)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

  vec K =y0.col(0)%th+(y0.col(1)%th%th)/2.0+(y0.col(2)%th%th%th)/6.0 +(y0.col(3)%th%th%th%th)/24.0;
  vec K1=y0.col(0)   +(y0.col(1)%th)       +(y0.col(2)%th%th)/2.0    +(y0.col(3)%th%th%th)/6.0;
  vec K2=y0.col(1)   +(y0.col(2)%th)       +(y0.col(3)%th%th)/2.0;
  vec val=exp(-0.5*log(2*3.141592653589793*K2)+(K-th%K1)-x0.col(8)-x0.col(9));

  p=(1.0/3.0) *(3*(y0.col(7)/6.0)%y0.col(5) - pow(y0.col(6)/2.0,2))/pow(y0.col(7)/6.0,2);
  q=(1.0/27.0)*(27*pow(y0.col(7)/6.0,2)%(y0.col(4)-Xt) - 9*(y0.col(7)/6.0)%(y0.col(6)/2.0)%y0.col(5) + 2*pow(y0.col(6)/2.0,3))/pow(y0.col(7)/6.0,3);
  chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  th=-(y0.col(6)/2.0)/(3*(y0.col(7)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

  K =y0.col(4)%th+(y0.col(5)%th%th)/2.0+(y0.col(6)%th%th%th)/6.0 +(y0.col(7)%th%th%th%th)/24.0;
  K1=y0.col(4)   +(y0.col(5)%th)       +(y0.col(6)%th%th)/2.0    +(y0.col(7)%th%th%th)/6.0;
  K2=y0.col(5)   +(y0.col(6)%th)       +(y0.col(7)%th%th)/2.0;
  val= log(val+exp(-0.5*log(2*3.141592653589793*K2)+(K-th%K1))%(1-exp(-x0.col(8)-x0.col(9))));

  List ret;
  ret["like"]  = val;
  ret["exc"]  = (1-exp(-x0.col(8)-x0.col(9)));
  ret["steps"]  = steps;
  return(ret);
}
  '
},{},{},{},{})}
 if(!factorize)
{
  switch(Dindex,
{
  Dpart='
  y0.col(0)= x0.col(0+4);
  y0.col(1)= x0.col(1+4)-1*y0.col(0)%x0.col(0+4);
  y0.col(2)= x0.col(2+4)-1*y0.col(0)%x0.col(1+4)-2*y0.col(1)%x0.col(0+4);
  y0.col(3)= x0.col(3+4)-1*y0.col(0)%x0.col(2+4)-3*y0.col(1)%x0.col(1+4)-3*y0.col(2)%x0.col(0+4);

  y0.col(4)= ((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
  y0.col(5)= ((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
  y0.col(6)= ((x0.col(2)-x0.col(6)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-2*y0.col(5)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));
  y0.col(7)= ((x0.col(3)-x0.col(7)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-1*y0.col(4)%((x0.col(2)-x0.col(6)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-3*y0.col(5)%((x0.col(1)-x0.col(5)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))))-3*y0.col(6)%((x0.col(0)-x0.col(4)%exp(-x0.col(8)-x0.col(9)))/(1-exp(-x0.col(8)-x0.col(9))));

  vec p=(1.0/3.0) *(3*(y0.col(3)/6.0)%y0.col(1) - pow(y0.col(2)/2.0,2))/pow(y0.col(3)/6.0,2);
  vec q=(1.0/27.0)*(27*pow(y0.col(3)/6.0,2)%(y0.col(0)-Xt) - 9*(y0.col(3)/6.0)%(y0.col(2)/2.0)%y0.col(1) + 2*pow(y0.col(2)/2.0,3))/pow(y0.col(3)/6.0,3);
  vec chk=pow(q,2)/4.0 + pow(p,3)/27.0;
  vec th=-(y0.col(2)/2.0)/(3*(y0.col(3)/6.0))+pow(-q/2.0+sqrt(chk),(1.0/3.0))-pow(q/2.0+sqrt(chk),(1.0/3.0));

  vec K =y0.col(0)%th+(y0.col(1)%th%th)/2.0+(y0.col(2)%th%th%th)/6.0 +(y0.col(3)%th%th%th%th)/24.0;
  vec K1=y0.col(0)   +(y0.col(1)%th)       +(y0.col(2)%th%th)/2.0    +(y0.col(3)%th%th%th)/6.0;
  vec K2=y0.col(1)   +(y0.col(2)%th)       +(y0.col(3)%th%th)/2.0;
  vec val=-0.5*log(2*3.141592653589793*K2)+(K-th%K1);

  List ret;
  ret["like"]  = val;
  ret["exc"]  = (1-exp(-x0.col(8)-x0.col(9)));
  ret["steps"]  = steps;
  return(ret);
}
  '
},{},{},{},{})}
}
if(DTR.order==6)
{

}
if(DTR.order==8)
{

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

   #if(state2)
   #{
   # DATA RESOLUTION -------------------------------------------------------------
   if((!homo.res)&(TR.order==4))
   {
    delt=cbind(diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,
               diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,diff(T.seq)/mesh,
               diff(T.seq)/mesh,diff(T.seq)/mesh)

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
                 c("1*(1+0*a.col(0))", "1*a.col(0)", "1*a.col(1)", ""                 , ""           , ""          ) ,
                 c("2*a.col(0)"      , "2*a.col(1)", "2*a.col(2)", "1*(1+0*a.col(0))", "1*a.col(0)" , "1*a.col(1)") ,
                 c("3*a.col(1)"      , "3*a.col(2)", "3*a.col(3)", "3*a.col(0)"       , "3*a.col(1)" , "3*a.col(2)") ,
                 c("4*a.col(2)"      , "4*a.col(3)", "4*(-(-5*a.col(3)%a.col(0)-10*a.col(2)%a.col(1)+20*a.col(2)%pow(a.col(0),2)+30*pow(a.col(1),2)%a.col(0)-60*a.col(1)%pow(a.col(0),3)+24*pow(a.col(0),5)))", "6*a.col(1)"       , "6*a.col(2)" , "6*a.col(3)") )
      MAT2=rbind(
                 c("1*(1+0*a.col(4))", "1*a.col(4)", "1*a.col(5)", ""                 , ""           , ""          ) ,
                 c("2*a.col(4)"      , "2*a.col(5)", "2*a.col(6)", "1*(1+0*a.col(4))", "1*a.col(4)" , "1*a.col(5)") ,
                 c("3*a.col(5)"      , "3*a.col(6)", "3*a.col(7)", "3*a.col(4)"       , "3*a.col(5)" , "3*a.col(6)") ,
                 c("4*a.col(6)"      , "4*a.col(7)", "4*(-(-5*a.col(7)%a.col(4)-10*a.col(6)%a.col(5)+20*a.col(6)%pow(a.col(4),2)+30*pow(a.col(5),2)%a.col(4)-60*a.col(5)%pow(a.col(4),3)+24*pow(a.col(4),5)))", "6*a.col(5)"       , "6*a.col(6)" , "6*a.col(7)") )
      MAT = rbind(MAT,MAT2)
   }
   if(TR.order==6)
   {
      MAT=rbind(
                 c("1*(1+0*a.col(0))", "1*a.col(0)", "1*a.col(1)", ""                 , ""           , ""          ) ,
                 c("2*a.col(0)"      , "2*a.col(1)", "2*a.col(2)", "1*(1+0*a.col(0))", "1*a.col(0)" , "1*a.col(1)") ,
                 c("3*a.col(1)"      , "3*a.col(2)", "3*a.col(3)", "3*a.col(0)"       , "3*a.col(1)" , "3*a.col(2)") ,
                 c("4*a.col(2)"      , "4*a.col(3)", "4*a.col(4)", "6*a.col(1)"       , "6*a.col(2)" , "6*a.col(3)") ,
                 c("5*a.col(3)"      , "5*a.col(4)", "5*a.col(5)", "10*a.col(2)"      , "10*a.col(3)", "10*a.col(4)"),
                 c("6*a.col(4)"      , "6*a.col(5)", "6*a.col(0)     ", "15*a.col(3)"      , "15*a.col(4)", "15*a.col(5)"))
   }
   if(TR.order==8)
   {
      MAT=rbind(
                 c("1*(1+0*a.col(0))", "1*a.col(0)", "1*a.col(1)", ""                 , ""           , ""          ) ,
                 c("2*a.col(0)"      , "2*a.col(1)", "2*a.col(2)", "1*(1+0*a.col(0))", "1*a.col(0)" , "1*a.col(1)") ,
                 c("3*a.col(1)"      , "3*a.col(2)", "3*a.col(3)", "3*a.col(0)"       , "3*a.col(1)" , "3*a.col(2)") ,
                 c("4*a.col(2)"      , "4*a.col(3)", "4*a.col(4)", "6*a.col(1)"       , "6*a.col(2)" , "6*a.col(3)") ,
                 c("5*a.col(3)"      , "5*a.col(4)", "5*a.col(5)", "10*a.col(2)"      , "10*a.col(3)", "10*a.col(4)"),
                 c("6*a.col(4)"      , "6*a.col(5)", "6*a.col(6)", "15*a.col(3)"      , "15*a.col(4)", "15*a.col(5)"),
                 c("7*a.col(5)"      , "7*a.col(6)", "7*a.col(7)", "21*a.col(4)"      , "21*a.col(5)", "21*a.col(6)"),
                 c("8*a.col(6)"      , "8*a.col(7)", "8*a.col(0)     ", "28*a.col(5)"      , "28*a.col(6)", "28*a.col(7)"))
 }

   # HOMOGENEITY -----------------------------------------------------------------
   namess2=c('G0','G1','G2','Q0','Q1','Q2','Lam0','Lam1','Jmu','Jsig','Jlam','Jalpha','Jbeta','Ja','Jb')
   func.list2=rep(0,length(namess2))
   obs=objects(pos=1)
   for(i in 1:length(namess2))
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

   # Jump Body Reference ---------------------------------------------------------
     jump.inhomogeneous=func.list.timehomo[-c(1:6)]
     mom.inhomogeneous =any(jump.inhomogeneous[-c(1,2)]==2)+1

   if(Jtype=='Mult')
   {
    jump.body =c(
                '+1*mm1',
                '+2*mm1+1*mm2',
                '+3*mm1+3*mm2+1*mm3',
                '+4*mm1+6*mm2+4*mm3+1*mm4',
                '+5*mm1+10*mm2+10*mm3+5*mm4+1*mm5',
                '+6*mm1+15*mm2+20*mm3+15*mm4+6*mm5+1*mm6',
                '+7*mm1+21*mm2+35*mm3+35*mm4+21*mm5+7*mm6+1*mm7',
                '+8*mm1+28*mm2+56*mm3+70*mm4+56*mm5+28*mm6+8*mm7+1*mm8')
    jump.body2=c(
                '+1*mm1',
                '+2*mm1+1*mm2',
                '+3*mm1+3*mm2+1*mm3',
                '+4*mm1+6*mm2+4*mm3+1*mm4',
                '+5*mm1+10*mm2+10*mm3+5*mm4+1*mm5',
                '+6*mm1+15*mm2+20*mm3+15*mm4+6*mm5+1*mm6',
                '+7*mm1+21*mm2+35*mm3+35*mm4+21*mm5+7*mm6+1*mm7',
                '+8*mm1+28*mm2+56*mm3+70*mm4+56*mm5+28*mm6+8*mm7+1*mm8')
    if(func.list2[7]==1)
    {
      for(i in 1:TR.order)
      {
         jump.body[i] = paste0('+a.col(',i-1,')',c('*','%')[jump.inhomogeneous[1]],'(',body(Lam0)[2],')',c('*','%')[mom.inhomogeneous],'(',jump.body[i],')')
      }
    }
    if(func.list2[8]==1)
    {
       for(i in 1:TR.order)
       {
          if(i !=4)
          {
          jump.body2[i] = paste0('+a.col(',i,')',c('*','%')[jump.inhomogeneous[2]],'(',body(Lam1)[2],')',c('*','%')[mom.inhomogeneous],'(',jump.body2[i],')')
          }
          if(i == 4){jump.body2[i] = paste0('+(-(-5*a.col(3)%a.col(0)-10*a.col(2)%a.col(1)+20*a.col(2)%pow(a.col(0),2)+30*pow(a.col(1),2)%a.col(0)-60*a.col(1)%pow(a.col(0),3)+24*pow(a.col(0),5)))',c('*','%')[jump.inhomogeneous[2]],'(',body(Lam1)[2],')',c('*','%')[mom.inhomogeneous],'(',jump.body2[i],')')
          }
       }
    }
   }
   if(Jtype=='Add')
   {
     if(mom.inhomogeneous==1)
     {
     jump.body =
              c("+1*mm1+0*a.col(0)"
              ,"+2*mm1*a.col(0)+1*mm2"
              ,"+3*mm1*a.col(1)+3*mm2*a.col(0)+1*mm3"
              ,"+4*mm1*a.col(2)+6*mm2*a.col(1)+4*mm3*a.col(0)+1*mm4"
              ,"+5*mm1*a.col(3)+10*mm2*a.col(2)+10*mm3*a.col(1)+5*mm4*a.col(0)+1*mm5"
              ,"+6*mm1*a.col(4)+15*mm2*a.col(3)+20*mm3*a.col(2)+15*mm4*a.col(1)+6*mm5*a.col(0)+1*mm6"
              ,"+7*mm1*a.col(5)+21*mm2*a.col(4)+35*mm3*a.col(3)+35*mm4*a.col(2)+21*mm5*a.col(1)+7*mm6*a.col(0)+1*mm7"
              ,"+8*mm1*a.col(6)+28*mm2*a.col(5)+56*mm3*a.col(4)+70*mm4*a.col(3)+56*mm5*a.col(2)+28*mm6*a.col(1)+8*mm7*a.col(0)+1*mm8")
     jump.body2 =
             c("+1*mm1*a.col(0)"
              ,"+2*mm1*a.col(1)+1*mm2*a.col(0)"
              ,"+3*mm1*a.col(2)+3*mm2*a.col(1)+1*mm3*a.col(0)"
              ,"+4*mm1*a.col(3)+6*mm2*a.col(2)+4*mm3*a.col(1)+1*mm4*a.col(0)"
              ,"+5*mm1*a.col(4)+10*mm2*a.col(3)+10*mm3*a.col(2)+5*mm4*a.col(1)+1*mm5*a.col(0)"
              ,"+6*mm1*a.col(5)+15*mm2*a.col(4)+20*mm3*a.col(3)+15*mm4*a.col(2)+6*mm5*a.col(1)+1*mm6*a.col(0)"
              ,"+7*mm1*a.col(6)+21*mm2*a.col(5)+35*mm3*a.col(4)+35*mm4*a.col(3)+21*mm5*a.col(2)+7*mm6*a.col(1)+1*mm7*a.col(0)"
              ,"+8*mm1*a.col(7)+28*mm2*a.col(6)+56*mm3*a.col(5)+70*mm4*a.col(4)+56*mm5*a.col(3)+28*mm6*a.col(2)+8*mm7*a.col(1)+1*mm8*a.col(0)")
     }
     if(mom.inhomogeneous==2)
     {
      jump.body =
             c("+1*mm1"
              ,"+2*mm1%a.col(0)+1*mm2"
              ,"+3*mm1%a.col(1)+3*mm2%a.col(0)+1*mm3"
              ,"+4*mm1%a.col(2)+6*mm2%a.col(1)+4*mm3%a.col(0)+1*mm4"
              ,"+5*mm1%a.col(3)+10*mm2%a.col(2)+10*mm3%a.col(1)+5*mm4%a.col(0)+1*mm5"
              ,"+6*mm1%a.col(4)+15*mm2%a.col(3)+20*mm3%a.col(2)+15*mm4%a.col(1)+6*mm5%a.col(0)+1*mm6"
              ,"+7*mm1%a.col(5)+21*mm2%a.col(4)+35*mm3%a.col(3)+35*mm4%a.col(2)+21*mm5%a.col(1)+7*mm6%a.col(0)+1*mm7"
              ,"+8*mm1%a.col(6)+28*mm2%a.col(5)+56*mm3%a.col(4)+70*mm4%a.col(3)+56*mm5%a.col(2)+28*mm6%a.col(1)+8*mm7%a.col(0)+1*mm8")
      jump.body2 =
             c("+1*mm1%a.col(0)"
              ,"+2*mm1%a.col(1)+1*mm2%a.col(0)"
              ,"+3*mm1%a.col(2)+3*mm2%a.col(1)+1*mm3%a.col(0)"
              ,"+4*mm1%a.col(3)+6*mm2%a.col(2)+4*mm3%a.col(1)+1*mm4%a.col(0)"
              ,"+5*mm1%a.col(4)+10*mm2%a.col(3)+10*mm3%a.col(2)+5*mm4%a.col(1)+1*mm5%a.col(0)"
              ,"+6*mm1%a.col(5)+15*mm2%a.col(4)+20*mm3%a.col(3)+15*mm4%a.col(2)+6*mm5%a.col(1)+1*mm6%a.col(0)"
              ,"+7*mm1%a.col(6)+21*mm2%a.col(5)+35*mm3%a.col(4)+35*mm4%a.col(3)+21*mm5%a.col(2)+7*mm6%a.col(1)+1*mm7%a.col(0)"
              ,"+8*mm1%a.col(7)+28*mm2%a.col(6)+56*mm3%a.col(5)+70*mm4%a.col(4)+56*mm5%a.col(3)+28*mm6%a.col(2)+8*mm7%a.col(1)+1*mm8%a.col(0)")
     }
     if(func.list2[7]==1)
     {
     for(i in 1:TR.order)
     {
        jump.body[i] = paste0('+(',body(Lam0)[2],')',c('*','%')[jump.inhomogeneous[1]],'(',jump.body[i],')')
     }
     }
     if(func.list2[8]==1)
     {
      for(i in 1:TR.order)
      {
         jump.body2[i] = paste0('+(',body(Lam1)[2],')',c('*','%')[jump.inhomogeneous[2]],'(',jump.body2[i],')')
      }
    }
   }



    prem =rep('   ',1)

    if(Jdist=='Exponential')
    {
        for(i in 1:TR.order)
        {
          prem =paste0(prem,c('double','vec')[jump.inhomogeneous[5]],' mm',i,'=',factorial(i),'*pow(',body(Jlam)[2],',',i,');\n    ')
        }
    }
    if(Jdist=='Normal')
    {
         prem =paste0(prem,c('double','vec')[jump.inhomogeneous[3]],' mu  =',body(Jmu)[2],';\n    ')
         prem =paste0(prem,c('double','vec')[jump.inhomogeneous[4]],' sig =',body(Jsig)[2],';\n    ')

           norm.prem =
         c(' mm1 = mu;'
        	,' mm2 = pow(mu,2)+ pow(sig,2);'
        	,' mm3 = pow(mu,3)+ 3* pow(mu,1)*pow(sig,2);'
        	,' mm4 = pow(mu,4)+ 6* pow(mu,2)*pow(sig,2)+3*pow(sig,4);'
        	,' mm5 = pow(mu,5)+ 10*pow(mu,3)*pow(sig,2)+15*pow(mu,1)*pow(sig,4);'
        	,' mm6 = pow(mu,6)+ 15*pow(mu,4)*pow(sig,2)+45*pow(mu,2)*pow(sig,4)+15*pow(sig,6);'
        	,' mm7 = pow(mu,7)+ 21*pow(mu,5)*pow(sig,2)+105*pow(mu,3)*pow(sig,4)+105*mu*pow(sig,6);'
        	,' mm8 = pow(mu,8)+ 28*pow(mu,6)*pow(sig,2)+210*pow(mu,4)*pow(sig,4)+420*pow(mu,2)*pow(sig,6)+105*pow(sig,8);')

         if((jump.inhomogeneous[3]==2)&&(jump.inhomogeneous[4]==2))
         {
           norm.prem =
         c(' mm1 = mu;'
        	,' mm2 = pow(mu,2)+ pow(sig,2);'
        	,' mm3 = pow(mu,3)+ 3* pow(mu,1)%pow(sig,2);'
        	,' mm4 = pow(mu,4)+ 6* pow(mu,2)%pow(sig,2)+3*pow(sig,4);'
        	,' mm5 = pow(mu,5)+ 10*pow(mu,3)%pow(sig,2)+15*pow(mu,1)%pow(sig,4);'
        	,' mm6 = pow(mu,6)+ 15*pow(mu,4)%pow(sig,2)+45*pow(mu,2)%pow(sig,4)+15*pow(sig,6);'
        	,' mm7 = pow(mu,7)+ 21*pow(mu,5)%pow(sig,2)+105*pow(mu,3)%pow(sig,4)+105*mu%pow(sig,6);'
        	,' mm8 = pow(mu,8)+ 28*pow(mu,6)%pow(sig,2)+210*pow(mu,4)%pow(sig,4)+420*pow(mu,2)%pow(sig,6)+105*pow(sig,8);')
         }

        for(i in 1:TR.order)
        {
          prem = paste0(prem,c('double','vec')[max(jump.inhomogeneous[3],jump.inhomogeneous[4])],norm.prem[i],'\n    ')
        }
    }

    if(Jdist=='Gamma')
    {
         ################################################################################################## Might be an error here  (%'s)
         prem =paste0(prem,c('double','vec')[jump.inhomogeneous[6]],' alphaa  =',body(Jalpha)[2],';\n    ')
         prem =paste0(prem,c('double','vec')[jump.inhomogeneous[7]],' betaa   =',body(Jbeta)[2],';\n    ')
         gam.prem =
         c(' mm1 = alphaa*betaa;'
        	,' mm2 = alphaa*(alphaa+1)*pow(betaa,2);'
        	,' mm3 = alphaa*(alphaa+1)*(alphaa+2)*pow(betaa,3);'
        	,' mm4 = alphaa*(alphaa+1)*(alphaa+2)*(alphaa+3)*pow(betaa,4);'
        	,' mm5 = 0;'
        	,' mm6 = 0;'
        	,' mm7 = 0;'
        	,' mm8 = 0;')
          if((jump.inhomogeneous[6]==2)&&(jump.inhomogeneous[7]==2))
         {
              gam.prem =
             c(' mm1 = alphaa%betaa;'
            	,' mm2 = alphaa%(alphaa+1)%pow(betaa,2);'
            	,' mm3 = alphaa%(alphaa+1)%(alphaa+2)%pow(betaa,3);'
            	,' mm4 = alphaa%(alphaa+1)%(alphaa+2)%(alphaa+3)%pow(betaa,4);'
            	,' mm5 = 0;'
            	,' mm6 = 0;'
            	,' mm7 = 0;'
            	,' mm8 = 0;')
         }

        for(i in 1:TR.order)
        {
          prem = paste0(prem,c('double','vec')[max(jump.inhomogeneous[6],jump.inhomogeneous[7])],gam.prem[i],'\n    ')
        }
    }

    if(Jdist=='Laplace')
    {
         ################################################################################################## Might be an error here (%'s)
         prem =paste0(prem,c('double','vec')[jump.inhomogeneous[8]],' aa  =',body(Ja)[2],';\n    ')
         prem =paste0(prem,c('double','vec')[jump.inhomogeneous[9]],' bb   =',body(Jb)[2],';\n    ')
         lap.prem =
         c(' mm1 = 0.5*(+2*aa*bb)                                                     ;'
        	,' mm2 = 0.5*(+2*pow(aa,2)+4*pow(bb,2))                                     ;'
        	,' mm3 = 0.5*(+2*pow(aa,3)+12*aa*pow(bb,2))                                  ;'
        	,' mm4 = 0.5*(+2*pow(aa,4)+24*pow(aa,2)*pow(bb,2)+48*pow(bb,4))               ;'
        	,' mm5 = 0.5*(+2*aa%bb)                                                     ;'
        	,' mm6 = 0.5*(+2*pow(aa,2)+4*pow(bb,2))                                     ;'
        	,' mm7 = 0.5*(+2*pow(aa,3)+12*aa%pow(bb,2))                                  ;'
        	,' mm8 = 0.5*(+2*pow(aa,4)+24*pow(aa,2)%pow(bb,2)+48*pow(bb,4))               ;')

         if((jump.inhomogeneous[8]==2)&&(jump.inhomogeneous[9]==2))
         {
              lap.prem =
               c(' mm1 = 0.5*(+2*aa%b)                                                     ;'
              	,' mm2 = 0.5*(+2*pow(aa,2)+4*pow(bb,2))                                     ;'
              	,' mm3 = 0.5*(+2*pow(aa,3)+12*aa%pow(bb,2))                                  ;'
              	,' mm4 = 0.5*(+2*pow(aa,4)+24*pow(aa,2)%pow(bb,2)+48*pow(bb,4))               ;'
              	,' mm5 = 0                                                                ;'
              	,' mm6 = 0                                                                ;'
              	,' mm7 = 0                                                                ;'
              	,' mm8 = 0                                                                ;')
         }

        for(i in 1:TR.order)
        {
          prem = paste0(prem,c('double','vec')[max(jump.inhomogeneous[8],jump.inhomogeneous[9])],lap.prem[i],'\n    ')
        }
    }

   #BUILD ODE --------------------------------------------------------------------
   dims=rep('(',TR.order*2)
   for(i in 1:TR.order)
   {
    for(j in which(func.list2[1:6]==1))
    {
      if(MAT[i,j]!='')
      {
        dims[i]=paste0(dims[i],'+(',body(namess2[j])[2],')',c('*','%')[func.list.timehomo[j]],'(',MAT[i,j],')')
        dims[i+TR.order]=paste0(dims[i+TR.order],'+(',body(namess2[j])[2],')',c('*','%')[func.list.timehomo[j]],'(',MAT[i+TR.order,j],')')

      }
    }
    if(func.list2[7]==1)
    {
      dims[i]=paste0(dims[i],jump.body[i])
    }
    if(func.list2[8]==1)
    {
      dims[i]=paste0(dims[i],jump.body2[i])
    }
    dims[i]=paste0(dims[i],')')
    dims[i+TR.order]=paste0(dims[i+TR.order],')')
   }

   if(any(dims=='()')){dims[which(dims=='()')]='0'}


   for(i in 1:(2*TR.order))
   {
      dims[i]=paste0(paste0(paste0('   atemp.col(',i-1,')='),dims[i]),';')
   }

   if(func.list2[7]==1)
   {
     jode1 = paste0('   atemp.col(',8,')=(',body(namess2[7])[2],')',c('*','%')[func.list.timehomo[7]],'(1+0*a.col(0));')
   }else
   {
     jode1 =  paste0('   atemp.col(',8,')=(0*a.col(0));')
   }
   if(func.list2[8]==1)
   {
     jode2 = paste0('   atemp.col(',9,')=(',body(namess2[8])[2],')',c('*','%')[func.list.timehomo[8]],'a.col(0);')
   }else
   {
     jode2 =  paste0('   atemp.col(',9,')=(0*a.col(0));')
   }


   # WRIGHT AND SOURCE -----------------------------------------------------------


   if(TR.order==4)
   {
     txt.full=paste(fpart,'\n',prem,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],'\n',dims[6],'\n',dims[7],'\n',dims[8],'\n',jode1,'\n',jode2,ODEpart,Dpart)
   }
   if(TR.order==6)
   {
     txt.full=paste(fpart,'\n',prem,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],'\n',dims[6],ODEpart,Dpart)
   }
   if(TR.order==8)
   {
     txt.full=paste(fpart,'\n',prem,'\n',dims[1],'\n',dims[2],'\n',dims[3],'\n',dims[4],'\n',dims[5],'\n',dims[6],'\n',dims[7],'\n',dims[8],ODEpart,Dpart)
   }
   type.sol ="          Jump Generalized Quadratic Diffusion (JGQD) "

   #library(Rcpp)
   #library(RcppArmadillo)
   if(wrt){write(txt.full,file='JGQD.mcmc.cpp')}
   sourceCpp(code=txt.full)
   #}
   #==============================================================================
   #                           Prior distribution Module
   #==============================================================================
   if(sum(objects(pos=1)=='priors')==1)
   {
     pp=function(theta){}
     body(pp)=parse(text =body(priors)[2])
     prior.list=paste0('d(theta)',':',paste0(body(priors)[2]))
   }
   if(sum(objects(pos=1)=='priors')==0)
   {
     prior.list=paste0('d(theta)',':',' No priors given.')
     pp=function(theta){1}
   }
   if(length(pp(theta))!=1){stop("Prior density must return only a single value!")}

   #==============================================================================
   #                           Interface Module
   #==============================================================================

    namess4=c('G0','G1','G2','Q0','Q1','Q2','Lam0','Lam1','Jmu','Jsig','Jlam','Jalpha','Jbeta','Ja','Jb')
    trim <- function (x) gsub("([[:space:]])", "", x)
    function.list=objects(pos=1)
    checklist=rep(0,length(namess4))
    for(i in 1:length(namess4))
    {
      checklist[i]= sum(function.list== namess4[i])
      if(sum(function.list==namess4[i]))
      {
        namess4[i]=paste0(namess4[i],' : ',trim(body(namess4[i])[2]))
      }
    }
    namess4=matrix(namess4,length(namess4),1)

    dinfo = c('Density approx. : ',

              'Trunc. Order    : ',
              'Dens.  Order    : ')
    dinfo[1] =paste0(dinfo[1],Dtype)

    dinfo[2] =paste0(dinfo[2],Trunc[1])
    dinfo[3] =paste0(dinfo[3],Trunc[2])

    buffer0=c('================================================================')
    buffer1=c('----------------------------------------------------------------')
    buffer2=c('................................................................')
    buffer3=c('...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ')
    buffer4=c('_____________________ Drift Coefficients _______________________')
    buffer5=c('___________________ Diffusion Coefficients _____________________')
    buffer6=c('__________________ Distribution Approximant ____________________')
    buffer7=c('_______________________ Jump Components ________________________')
    buffer8=c('......................... Intensity ............................')
    buffer9=c('........................... Jumps ..............................')
    type.sol ="           Jump Generalized Quadratic Diffusion (JGQD) "



   if(Jdist=='Normal')     {checklist[12:16-1]=0}
   if(Jdist=='Exponential'){checklist[c(c(10,11,13,14,15,16)-1)]=0}
   if(Jdist=='Gamma')      {checklist[c(10:12,15,16)-1]=0}
   if(Jdist=='Laplace')      {checklist[c(10:14)-1]=0}
    Info=c(buffer0,type.sol,buffer0,buffer4,namess4[1:3],buffer5,namess4[4:6],buffer7,namess4[7:8],buffer9,Jdist,namess4[8+which(checklist[-(1:8)]==1)],buffer6,dinfo)
    Info=data.frame(matrix(Info,length(Info),1))
    colnames(Info)=''
    if(print.output)
    {print(Info,row.names = FALSE,right=F)}

    ############################################################################
    ############################################################################
    ############################################################################

    tme=Sys.time()


    par.matrix=matrix(0,length(theta),updates)
    ll=rep(0,updates)
    acc=ll
    kk=0
    par.matrix[,1]=theta
    prop.matrix =par.matrix
    rs=solver(X[-nnn],X[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],P,alpha,lower,upper,TR.order)
    lold=sum(rs$like)
    ll[1]=lold

    excess = matrix(0,2,updates)
    term.discont  = rep(0,nnn-1)
    oldexc = term.discont
     pb <- txtProgressBar(1,updates,1,style = 1,width = 65)
    failed.chain=F
    retries = 0
    retry.count   = 0
    retry.indexes = c()
    max.retries = 0
    success = TRUE
    stps= rep(0,updates)
    i = 2

    while(i<=updates)
    {
        theta.temp=theta

        theta=theta+rnorm(length(theta),sd=sds)
        prop.matrix[,i] = theta
        rs=solver(X[-nnn],X[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],P,alpha,lower,upper,TR.order)
        lnew=sum(rs$like)
        stps[i]=(rs$steps)
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
            rs=solver(X[-nnn],X[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],P,alpha,lower,upper,TR.order)
            lnew=sum(rs$like)
            stps[i]=(rs$steps)
            rat=min(exp(lnew-lold)*pp(theta)/pp(theta.temp),1)
            if(is.na(rat)){retry.count=retry.count+1}
          }
        }
        u=runif(1)
        is.true =(rat>u)
        is.false=!is.true
        excess[1,i] = mean(rs$exc)*is.true+excess[1,i-1]*is.false
        oldexc= rs$exc*is.true+oldexc*is.false
        term.discont =term.discont+oldexc
        theta=theta*is.true+theta.temp*is.false
        lold=lnew*is.true +lold*is.false

        par.matrix[,i]=theta
        ll[i]=lold
        kk=kk+is.true
        acc[i]=is.true
        if(any(is.na(theta))){print('Fail'); ;failed.chain=TRUE;break;}
         if(max.retries>5000){print('Fail: Failed evaluation limit exceeded!');failed.chain=TRUE;break;}
        setTxtProgressBar(pb, i)
        i=i+1
    }
    close(pb)
    acc = cumsum(acc)/(1:updates)
    tme.eval = function(start_time)
    {
      start_time = as.POSIXct(start_time)
      dt = difftime(Sys.time(), start_time, units="secs")
      format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
     }
    tme=tme.eval(tme)

    term.discont =term.discont/updates

    theta =apply(par.matrix[,-c(1:burns)],1,mean)
    meanD=mean(-2*ll[-c(1:burns)])
    pd=meanD-(-2*sum(solver(X[-nnn],X[-1],c(0,theta),mesh,delt,nnn-1,T.seq[-nnn],P,alpha,lower,upper,TR.order)$like))
    DIC=meanD+pd
    actual.p=length(theta)

    model.inf=list(elapsed.time=tme,time.homogeneous=c('Yes','No')[2-homo],p=actual.p,DIC=DIC,pd=pd,N=length(X)-length(excl)+1,Tag=Tag)
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
      d1=1:((nper)+2)
      d2=d1
      O=outer(d1,d2)
      test=O-((nper)+2)
      test[test<0]=100
      test=test[1:4,1:4]
      wh=which(test==min(test))
      wh
      d1=d1[col(test)[wh[1]]]
      d2=d2[row(test)[wh[1]]]
      par(mfrow=c(d1,d2))
      #plot(pmin(stps,50),type='s')
      cols=rainbow_hcl(nper, start = 10, end = 275,c=100,l=70)
      ylabs=paste0('theta[',1:nper,']')
      for(i in 1:nper)
      {
          plot(prop.matrix[i,],col='gray90',type='s',main=ylabs[i],xlab='Iteration',ylab='')
          lines(par.matrix[i,],col=cols[i],type='s')
          abline(v=burns,lty='dotdash')

      }
      plot(acc,type='l',ylim=c(0,1),col='darkblue',main='Accept. Rate',xlab='Iteration',ylab='%/100')
      abline(h=seq(0,1,1/10),lty='dotted')
      abline(v=burns,lty='dotdash')
      abline(h=0.4,lty='solid',col='red',lwd=1.2)
      abline(h=0.2,lty='solid',col='red',lwd=1.2)
      if(length(retry.indexes)>0)
      {
          axis(1,at=retry.indexes,labels =NA,tcl=-0.2,col='grey')
      }
      box()
      #windows()
      #par(mfrow=c(1,2))
      lines(excess[1,],col='gray')
      lines(acc,col='darkblue')
      #abline(h=seq(0,1,1/10),lty='dotted')
      plot(term.discont~time[-1],type='l',col='darkblue',ylim=c(0,max(0.5,round(term.discont,2))))
      abline(h=seq(0,1,1/10),lty='dotted')
      lines((abs(diff(X))/max(abs(diff(X)))*max(0.5,round(term.discont,2)))~time[-1],col='gray75',type='h')
      lines(term.discont~time[-1],col='darkblue')

    }
    ret=list(par.matrix=t(par.matrix),acceptance.rate=acc,elapsed.time=tme,model.info=model.inf,failed.chain=failed.chain,retry.indexes=retry.indexes,zero.jump=excess[1,],prop.matrix=t(prop.matrix))
    class(ret) = 'JGQD.mcmc'
    return(ret)
  }
