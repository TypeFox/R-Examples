globalVariables(c('mua','alph'))
JGQD.density=function(Xs=4,Xt=seq(5,8,1/10),s=0,t=5,delt=1/100,Jdist='Normal',Jtype='Add',Dtype='Saddlepoint',Trunc=c(8,4),factorize=FALSE,factor.type='Diffusion',beta,print.output=TRUE,eval.density=TRUE)
{
   Integ = TRUE
   P=100;alpha=0.1;lower=0;upper=50
   TR.order=Trunc[1]
   DTR.order=Trunc[2]
   Dtypes =c('Saddlepoint','Edgeworth','Normal.A','Gamma.A')

   check_for_model=function()
  {
    txt=''
    namess=c('G0','G1','G2','Q0','Q1','Q2','Lam0','Lam1','Lam2','Jmu','Jsig','Jlam','Jalpha','Jbeta','Ja','Jb')
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
  G0=function(t){5+2*sin(4*pi*t)}
  G1=function(t){-1}
  Q1=function(t){1}

  model=GQD.density(5,seq(0,10,1/10),0,5,1/100)
  GQD.plot(model)
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
  Q0=function(t){2+sin(5*pi*t)}
  model=GQD.density(5,seq(0,10,1/10),0,5,1/100)
  GQD.plot(model)
  --------------------------------------------------------------------------------
   '
     check=T
    }
    return(list(check=check,txt=txt))
  }
  check_for=check_for_model()
  if(check_for[[1]]){stop(check_for[[2]])}

  TR.order=Trunc[1]
  DTR.order=Trunc[2]
  #Dtypes =c('Saddlepoint','Normal','Beta')
  Dindex = which(Dtypes==Dtype)
  JDtypes =c('Normal','Exponential','Gamma','Laplace','Uniform')
  JDindex = which(JDtypes==Jdist)

  b1 = '\n==============================================================================\n'
  b2 = '==============================================================================\n'
  warn=c(
     'Input  1: Argument {lower} must be of length 1.\n'
    ,'Input  2: Argument {upper} must be of length 1.\n'
    ,'Input  3: {upper} must be > {lower}.\n'
    ,'Input  4: Range [lower,upper] must be strictly positive for Dtype Gamma or InvGamma.\n'
    ,'Input  5: Dtype cannot be Beta for observations not in (0,1).\n'
    ,'Input  6: P must be >= 10.\n'
    ,'Input  7: Dtype has to be one of Saddlepoint or Edgeworth.\n'
    ,'Input  8: Trunc[2] must be <= Trunc[1].\n'
    ,'Input  9: Large {delt} may result in poor approximations.\n'
    ,'Input 10: Argument {delt} must be < 1.\n'
    ,'Input 11: length(P)!=1.\n'
    ,'Input 12: length(delt)!=1.\n'
    ,'Input 13: length(alpha)!=1.\n'
    ,'Input 14: length(Xt) must be > 1.\n'
    ,'Input 15: Argument {Xs} must be of length 1.\n'
    ,'Input 16: Argument {Xt} must be of type vector!.\n'
    ,'Input 17: Argument {time} must be of type vector!.\n'
    ,'Input 18: Argument {alpha} must be of length 1.\n'
    ,'Input 19: Argument {Trunc} must be of length 2.\n'
    ,'Input 20: Cumulant truncation (Trunc[1]) must be either 4 or 8.\n'
    ,'Input 21: Density truncation (Trunc[2]) must be either 4 or 8.\n'
    ,'Input 22: Starting time {s} cannot be greater than {t}.\n'
    ,'Input 23: {Jdist} must be one of Normal, Exponential, Gamma or Laplace.\n'
    ,'Input 24: {Jtype} has to be of type Add or Mult.\n'
    ,'Input 25: {factorize} has to be TRUE or FALSE.\n'
    ,'Input 26: {factor.type} has to be Diffusion or Hawke.\n'
    ,'Input 27: Diffusion terms are 0. Use {factorize=TRUE} and {factor.type="Hawke"}.\n'
    ,'Input 28: Diffusion terms are 0. Density approximant cab perform poorly without \nfactorization. Try {factorize=TRUE}.\n'

  )
  warntrue=rep(FALSE,40)

  if(length(Xt)<2){warntrue[14] =TRUE}#{stop(paste0(b1,warn[14],b2))}
  if(length(Xs)!=1){warntrue[15] =TRUE}#{stop(paste0(b1,warn[15],b2))}
  if(!is.vector(Xt)){warntrue[16] =TRUE}#{stop(paste0(b1,warn[16],b2))}
  if(sum(Dindex)==0){warntrue[7] =TRUE}#{stop(paste0(b1,warn[7],b2))}
  if(length(lower)>1){warntrue[1] =TRUE}#{stop(paste0(b1,warn[1],b2))}
  if(length(upper)>1){warntrue[2] =TRUE}#{stop(paste0(b1,warn[2],b2))}
  if(upper<=lower)   {warntrue[3] =TRUE}# {stop(paste0(b1,warn[3],b2))}
  #if((Dindex==3)|(Dindex==4)){if(lower[1]<=0){warntrue[4] =TRUE}}#{stop(paste0(b1,warn[4],b2))}}
  if(Dindex==5){if(any(Xt<=0)|any(Xt>=1)){warntrue[5] =TRUE}}#{stop(paste0(b1,warn[5],b2))}}
  if(P<10){warntrue[6] =TRUE}#{stop(paste0(b1,warn[6],b2))}
  if(length(P)!=1){warntrue[11] =TRUE}#{stop(paste0(b1,warn[11],b2))}
  if(length(alpha)!=1){warntrue[18] =TRUE}#{stop(paste0(b1,warn[18],b2))}
  if(Trunc[2]>Trunc[1]){warntrue[8] =TRUE}#{stop(paste0(b1,warn[8],b2))}
  if(length(Trunc)!=2){warntrue[19] =TRUE}#{stop(paste0(b1,warn[19],b2))}
  if(sum(TR.order==c(4,8))!=1){warntrue[20] =TRUE}#{stop(paste0(b1,warn[20],b2))}
  if(sum(DTR.order==c(4,8))!=1){warntrue[21] =TRUE}#{stop(paste0(b1,warn[21],b2))}
  #if(delt>=1/5){stop(paste0(b1,warn[9],b2))}
  if(delt>=1){warntrue[10] =TRUE}#{stop(paste0(b1,warn[10],b2))}
  if(t<s){warntrue[22] =TRUE}#{stop(paste0(b1,warn[22],b2))}
  if((Jtype!='Add')&&(Jtype!='Mult'))                          {warntrue[24] =TRUE}
  if((factorize!=TRUE)&&(factorize!=FALSE))                    {warntrue[25] =TRUE}
  if((factor.type!='Diffusion')&&(factor.type!='Hawke'))       {warntrue[26] =TRUE}
  if(sum(JDindex)==0){stop(paste0(b1,warn[23],b2))}

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

   diffusion.function=function()
 {
   time =seq(s,t,delt)[-1]
   bod='0*x^0'
   coefs=c('Q0','Q1','Q2')
   obs=objects(pos=1)
   for(i in 1:length(coefs))
   {
    if(sum(obs==coefs[i])){bod=paste0(bod,'+(',body(coefs[i])[2],')*x^',i-1)}
   }
   difff=function(x,t){}
   body(difff)  = parse(text=bod)
   #print(body(difff))
   M = matrix(0,length(Xt),length(time))
   for(i in 1:length(time))
   {
      M[,i]=round(difff(Xt,time[i]),6)>0
   }
   if(sum(M)==length(time)*length(Xt)){return(1)}
   if(sum(M)!=0)                      {return(2)}
   if(sum(M)==0)                      {return(3)}
 }
 switch(diffusion.function(),
 {},
 {},
 {
   if((factorize)&&(factor.type!='Hawke')){stop(paste0(b1,warn[27],b2))}
   if(!factorize){warning(paste0(b1,warn[28],b2))}
 })

  pow=function(x,p)
  {
    x^p
  }
  prod=function(a,b){a*b}

   function.list=objects(pos=1)
   namess=c('G0','G1','G2','Q0','Q1','Q2','Lam0','Lam1','Lam2','Jmu','Jsig','Jlam','Jalpha','Jbeta','Ja','Jb')
   coef.terms2=paste0(namess,'(t)')
   checklist=rep(0,length(namess))
   for(i in 1:length(namess))
   {
        checklist[i]= sum(function.list== namess[i])
   }
   if(factorize)
   {
    if(sum(checklist[c(3,5,6)])==0){Dtype='Normal.A'}
   }

   if(Jdist=='Normal')     {checklist[12:16]=0}
   if(Jdist=='Exponential'){checklist[c(10,11,13,14,15,16)]=0}
   if(Jdist=='Gamma')      {checklist[c(10:12,15,16)]=0}
   if(Jdist=='Laplace')      {checklist[c(10:14)]=0}
   if(Jdist=='Uniform')      {checklist[c(10:14)]=0}
   coef.index=which(checklist==1)

    #==============================================================================
  #                           Interface Module
  #==============================================================================

  namess4=namess
  trim <- function (x) gsub("([[:space:]])", "", x)

  for(i in 1:length(namess4))
  {
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
  #if(Dindex!=1)
  #{
  #  dinfo[2] =paste0(dinfo[2],P)
  #  dinfo[3] =paste0(dinfo[3],alpha)
  #}
  dinfo[2] =paste0(dinfo[2],Trunc[1])
  dinfo[3] =paste0(dinfo[3],Trunc[2])

  buffer0=c('================================================================')
  buffer1=c('----------------------------------------------------------------')
  buffer2=c('................................................................')
  buffer3=c('...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ')
  buffer4=c('_____________________ Drift Coefficients _______________________')
  buffer5=c('___________________ Diffusion Coefficients _____________________')
  buffer6=c('__________________ Distribution Approximant ____________________')
  buffer7=c('_______________________ Jump Mechanism _________________________')
  buffer8=c('......................... Intensity ............................')
  buffer9=c('........................... Jumps ..............................')

  type.sol ="           Jump Generalized Quadratic Diffusion (JGQD) "

  Info=c(buffer0,type.sol,buffer0,buffer4,namess4[1:3],buffer5,namess4[4:6],buffer7,buffer8,namess4[7:9],buffer9,Jdist,namess4[9+which(checklist[-(1:9)]==1)],buffer6,dinfo)
  Info=data.frame(matrix(Info,length(Info),1))
  colnames(Info)=''

  if(print.output)
  {
    print(Info,row.names = FALSE,right=F)
  }

if(Integ)
{
  if(TR.order ==8)
 {
   MAT=rbind(
    c('1','1*x[1]','1*x[2]','','',''),
    c('2*x[1]','2*x[2]','2*x[3]','1','1*x[1]','1*x[2]'),
    c('3*x[2]','3*x[3]','3*x[4]','3*x[1]','3*x[2]','3*x[3]'),
    c('4*x[3]','4*x[4]','4*x[5]','6*x[2]','6*x[3]','6*x[4]'),
    c('5*x[4]','5*x[5]','5*x[6]','10*x[3]','10*x[4]','10*x[5]'),
    c('6*x[5]','6*x[6]','6*x[7]','15*x[4]','15*x[5]','15*x[6]'),
    c('7*x[6]','7*x[7]','7*x[8]','21*x[5]','21*x[6]','21*x[7]'),
    c('8*x[7]','8*x[8]','8*(-(-9*x[1]*x[8]-36*x[2]*x[7]+72*x[1]^2*x[7]-84*x[3]*x[6]+504*x[1]*x[2]*x[6]-504*x[1]^3*x[6]-126*x[4]*x[5]+1008*x[1]*x[3]*x[5]+756*
x[2]^2*x[5]-4536*x[1]^2*x[2]*x[5]+3024*x[1]^4*x[5]+630*x[1]*x[4]^2+2520*x[2]*x[3]*x[4]-7560*x[1]^2*x[3]*x[4]-11340*x[1]*x[2]^2*x[4]+30240*x[1]^3*x[2]*x[4]
-15120*x[1]^5*x[4]+560*x[3]^3-15120*x[1]*x[2]*x[3]^2+20160*x[1]^3*x[3]^2-7560*x[2]^3*x[3]+90720*x[1]^2*x[2]^2*x[3]-151200*x[1]^4*x[2]*x[3]+60480*
x[1]^6*x[3]+22680*x[1]*x[2]^4-151200*x[1]^3*x[2]^3+272160*x[1]^5*x[2]^2-181440*x[1]^7*x[2]+40320*x[1]^9))','28*x[6]','28*x[7]','28*x[8]'))

MAT2=rbind(
    c('1','1*x[8+1]','1*x[8+2]','','',''),
    c('2*x[8+1]','2*x[8+2]','2*x[8+3]','1','1*x[8+1]','1*x[8+2]'),
    c('3*x[8+2]','3*x[8+3]','3*x[8+4]','3*x[8+1]','3*x[8+2]','3*x[8+3]'),
    c('4*x[8+3]','4*x[8+4]','4*x[8+5]','6*x[8+2]','6*x[8+3]','6*x[8+4]'),
    c('5*x[8+4]','5*x[8+5]','5*x[8+6]','10*x[8+3]','10*x[8+4]','10*x[8+5]'),
    c('6*x[8+5]','6*x[8+6]','6*x[8+7]','15*x[8+4]','15*x[8+5]','15*x[8+6]'),
    c('7*x[8+6]','7*x[8+7]','7*x[8+8]','21*x[8+5]','21*x[8+6]','21*x[8+7]'),
    c('8*x[8+7]','8*x[8+8]','8*(-(-9*x[8+1]*x[8+8]-36*x[8+2]*x[8+7]+72*x[8+1]^2*x[8+7]-84*x[8+3]*x[8+6]+504*x[8+1]*x[8+2]*x[8+6]-504*x[8+1]^3*x[8+6]-126*x[8+4]*x[8+5]+1008*x[8+1]*x[8+3]*x[8+5]+756*
x[8+2]^2*x[8+5]-4536*x[8+1]^2*x[8+2]*x[8+5]+3024*x[8+1]^4*x[8+5]+630*x[8+1]*x[8+4]^2+2520*x[8+2]*x[8+3]*x[8+4]-7560*x[8+1]^2*x[8+3]*x[8+4]-11340*x[8+1]*x[8+2]^2*x[8+4]+30240*x[8+1]^3*x[8+2]*x[8+4]
-15120*x[8+1]^5*x[8+4]+560*x[8+3]^3-15120*x[8+1]*x[8+2]*x[8+3]^2+20160*x[8+1]^3*x[8+3]^2-7560*x[8+2]^3*x[8+3]+90720*x[8+1]^2*x[8+2]^2*x[8+3]-151200*x[8+1]^4*x[8+2]*x[8+3]+60480*
x[8+1]^6*x[8+3]+22680*x[8+1]*x[8+2]^4-151200*x[8+1]^3*x[8+2]^3+272160*x[8+1]^5*x[8+2]^2-181440*x[8+1]^7*x[8+2]+40320*x[8+1]^9))','28*x[8+6]','28*x[8+7]','28*x[8+8]'))

   dims=rep('',8)
   dims2=rep('',8)
   for(i in which(checklist[1:6]==1))
   {
     for(j in 1:8)
     {
      if(MAT[j,i]!='')
      {
                dims[j] = paste0(dims[j],'+(',body(namess[i])[2],')*',MAT[j,i])
                dims2[j] = paste0(dims2[j],'+(',body(namess[i])[2],')*',MAT2[j,i])
      }
     }
   }

   if(Jtype=='Add')
   {
        JA0 =
        c("+1*mm1"
        ,"+2*mm1*x[1]+1*mm2"
        ,"+3*mm1*x[2]+3*mm2*x[1]+1*mm3"
        ,"+4*mm1*x[3]+6*mm2*x[2]+4*mm3*x[1]+1*mm4"
        ,"+5*mm1*x[4]+10*mm2*x[3]+10*mm3*x[2]+5*mm4*x[1]+1*mm5"
        ,"+6*mm1*x[5]+15*mm2*x[4]+20*mm3*x[3]+15*mm4*x[2]+6*mm5*x[1]+1*mm6"
        ,"+7*mm1*x[6]+21*mm2*x[5]+35*mm3*x[4]+35*mm4*x[3]+21*mm5*x[2]+7*mm6*x[1]+1*mm7"
        ,"+8*mm1*x[7]+28*mm2*x[6]+56*mm3*x[5]+70*mm4*x[4]+56*mm5*x[3]+28*mm6*x[2]+8*mm7*x[1]+1*mm8")

        JA1 =
        c("+1*mm1*x[1]"
         ,"+2*mm1*x[2]+1*mm2*x[1]"
         ,"+3*mm1*x[3]+3*mm2*x[2]+1*mm3*x[1]"
         ,"+4*mm1*x[4]+6*mm2*x[3]+4*mm3*x[2]+1*mm4*x[1]"
         ,"+5*mm1*x[5]+10*mm2*x[4]+10*mm3*x[3]+5*mm4*x[2]+1*mm5*x[1]"
         ,"+6*mm1*x[6]+15*mm2*x[5]+20*mm3*x[4]+15*mm4*x[3]+6*mm5*x[2]+1*mm6*x[1]"
         ,"+7*mm1*x[7]+21*mm2*x[6]+35*mm3*x[5]+35*mm4*x[4]+21*mm5*x[3]+7*mm6*x[2]+1*mm7*x[1]"
         ,"+8*mm1*x[8]+28*mm2*x[7]+56*mm3*x[6]+70*mm4*x[5]+56*mm5*x[4]+28*mm6*x[3]+8*mm7*x[2]+1*mm8*x[1]")

                 JA2 =
        c("+1*mm1*x[2]"
         ,"+2*mm1*x[3]+1*mm2*x[2]"
         ,"+3*mm1*x[4]+3*mm2*x[3]+1*mm3*x[2]"
         ,"+4*mm1*x[5]+6*mm2*x[4]+4*mm3*x[3]+1*mm4*x[2]"
         ,"+5*mm1*x[6]+10*mm2*x[5]+10*mm3*x[4]+5*mm4*x[3]+1*mm5*x[2]"
         ,"+6*mm1*x[7]+15*mm2*x[6]+20*mm3*x[5]+15*mm4*x[4]+6*mm5*x[3]+1*mm6*x[2]"
         ,"+7*mm1*x[8]+21*mm2*x[7]+35*mm3*x[6]+35*mm4*x[5]+21*mm5*x[4]+7*mm6*x[3]+1*mm7*x[2]"
         ,"+8*mm1*((-(-9*x[8+1]*x[8+8]-36*x[8+2]*x[8+7]+72*x[8+1]^2*x[8+7]-84*x[8+3]*x[8+6]+504*x[8+1]*x[8+2]*x[8+6]-504*x[8+1]^3*x[8+6]-126*x[8+4]*x[8+5]+1008*x[8+1]*x[8+3]*x[8+5]+756*
x[8+2]^2*x[8+5]-4536*x[8+1]^2*x[8+2]*x[8+5]+3024*x[8+1]^4*x[8+5]+630*x[8+1]*x[8+4]^2+2520*x[8+2]*x[8+3]*x[8+4]-7560*x[8+1]^2*x[8+3]*x[8+4]-11340*x[8+1]*x[8+2]^2*x[8+4]+30240*x[8+1]^3*x[8+2]*x[8+4]
-15120*x[8+1]^5*x[8+4]+560*x[8+3]^3-15120*x[8+1]*x[8+2]*x[8+3]^2+20160*x[8+1]^3*x[8+3]^2-7560*x[8+2]^3*x[8+3]+90720*x[8+1]^2*x[8+2]^2*x[8+3]-151200*x[8+1]^4*x[8+2]*x[8+3]+60480*
x[8+1]^6*x[8+3]+22680*x[8+1]*x[8+2]^4-151200*x[8+1]^3*x[8+2]^3+272160*x[8+1]^5*x[8+2]^2-181440*x[8+1]^7*x[8+2]+40320*x[8+1]^9)))  +28*mm2*x[8]+56*mm3*x[7]+70*mm4*x[6]+56*mm5*x[5]+28*mm6*x[4]+8*mm7*x[3]+1*mm8*x[2]")

        if(checklist[7]==1)
        {
          for(i in 1:8)
          {
             dims[i] = paste0(dims[i],'+(',body(namess[7])[2],')*(',JA0[i],')')
          }
        }
        if(checklist[8]==1)
        {
           for(i in 1:8)
          {
             dims[i] = paste0(dims[i],'+(',body(namess[8])[2],')*(',JA1[i],')')
          }
        }
        if(checklist[9]==1)
        {
           for(i in 1:8)
          {
             dims[i] = paste0(dims[i],'+(',body(namess[9])[2],')*(',JA2[i],')')
          }
        }

   }
   #c('G0','G1','G2','Q0','Q1','Q2','Lam0','Lam1','Lam2','Jmu','Jsig','Jlam','Jalpha','Jbeta')
   if(checklist[7]==0){Lam1=function(t){0}}
   if(checklist[8]==0){Lam1=function(t){0}}
   if(checklist[9]==0){Lam2=function(t){0}}
   if(checklist[10]==0){Jmu=function(t){0}}
   if(checklist[11]==0){Jsig=function(t){0}}
   if(checklist[12]==0){Jlam=function(t){0}}
   if(checklist[13]==0){Jalpha=function(t){0}}
   if(checklist[14]==0){Jbeta=function(t){0}}
   if(checklist[15]==0){Ja=function(t){0}}
   if(checklist[16]==0){Jb=function(t){0}}

   odekernel=paste0('c(',dims[1],'\n,',dims[2],'\n,',dims[3],'\n,',dims[4],'\n,',dims[5],'\n,',dims[6],'\n,',dims[7],'\n,',dims[8],'\n,',
   dims2[1],'\n,',dims2[2],'\n,',dims2[3],'\n,',dims2[4],'\n,',dims2[5],'\n,',dims2[6],'\n,',dims2[7],'\n,',dims2[8],'\n,',
   body(namess[7])[2],'\n,',paste('(',body(namess[8])[2],')*x[9]'),'\n,',paste('(',body(namess[9])[2],')*x[10]'),')')
   }
     if(TR.order ==4)
 {
     MAT=rbind(
    c('1','1*x[1]','1*x[2]','','',''),
    c('2*x[1]','2*x[2]','2*x[3]','1','1*x[1]','1*x[2]'),
    c('3*x[2]','3*x[3]','3*x[4]','3*x[1]','3*x[2]','3*x[3]'),
    c('4*x[3]','4*x[4]','4*x[5]','6*x[2]','6*x[3]','6*x[4]'))

    MAT2=rbind(
    c('1','1*x[4+1]','1*x[4+2]','','',''),
    c('2*x[4+1]','2*x[4+2]','2*x[4+3]','1','1*x[4+1]','1*x[4+2]'),
    c('3*x[4+2]','3*x[4+3]','3*x[4+4]','3*x[4+1]','3*x[4+2]','3*x[4+3]'),
    c('4*x[4+3]','4*x[4+4]','4*x[4+5]','6*x[4+2]','6*x[4+3]','6*x[4+4]'))

   dims=rep('',4)
   dims2=rep('',4)
   for(i in which(checklist[1:6]==1))
   {
     for(j in 1:4)
     {
      if(MAT[j,i]!='')
      {
                dims[j] = paste0(dims[j],'+(',body(namess[i])[2],')*',MAT[j,i])
                dims2[j] = paste0(dims2[j],'+(',body(namess[i])[2],')*',MAT2[j,i])
      }
     }
   }

   if(Jtype=='Add')
   {
        JA0 =
        c("+1*mm1"
        ,"+2*mm1*x[1]+1*mm2"
        ,"+3*mm1*x[2]+3*mm2*x[1]+1*mm3"
        ,"+4*mm1*x[3]+6*mm2*x[2]+4*mm3*x[1]+1*mm4")

        JA1 =
        c("+1*mm1*x[1]"
         ,"+2*mm1*x[2]+1*mm2*x[1]"
         ,"+3*mm1*x[3]+3*mm2*x[2]+1*mm3*x[1]"
         ,"+4*mm1*x[4]+6*mm2*x[3]+4*mm3*x[2]+1*mm4*x[1]")

                 JA2 =
        c("+1*mm1*x[2]"
         ,"+2*mm1*x[3]+1*mm2*x[2]"
         ,"+3*mm1*x[4]+3*mm2*x[3]+1*mm3*x[2]"
         ,"+4*mm1*x[5]+6*mm2*x[4]+4*mm3*x[3]+1*mm4*x[2]")

   }
    if((!missing(beta))&&(length(beta)==2))
   {
        a=beta[1]
        b=beta[2]
        JA0 =
        c("b*x[1]*mm1+a*mm1"
        ,"b^2*x[2]*mm2+2*a*b*x[1]*mm2+a^2*mm2+2*b*x[2]*mm1+2*a*x[1]*mm1"
        ,"b^3*x[3]*mm3+3*a*b^2*x[2]*mm3+3*a^2*b*x[1]*mm3+a^3*mm3+3*b^2*x[3]*mm2+6*a*b*x[2]*mm2+3*a^2*x[1]*mm2+3*b*x[3]*mm1+3*a*x[2]*mm1"
        ,"b^4*x[4]*mm4+4*a*b^3*x[3]*mm4+6*a^2*b^2*x[2]*mm4+4*a^3*b*x[1]*mm4+a^4*mm4+4*b^3*x[4]*mm3+12*a*b^2*x[3]*mm3+12*a^2*b*x[2]*mm3+4*a^3*x[1]*mm3+6*b^2*x[4]*mm2+12*a*b*x[3]*mm2+6*a^2*x[2]*mm2+4*b*x[4]*mm1+4*a*x[3]*mm1")
        print('Invoked special jump structure.')
   }
  if(checklist[7]==1)
  {
    for(i in 1:4)
    {
       dims[i] = paste0(dims[i],'+(',body(namess[7])[2],')*(',JA0[i],')')
    }
  }
  if(checklist[8]==1)
  {
     for(i in 1:4)
    {
       dims[i] = paste0(dims[i],'+(',body(namess[8])[2],')*(',JA1[i],')')
    }
  }
  if(checklist[9]==1)
  {
     for(i in 1:4)
    {
       dims[i] = paste0(dims[i],'+(',body(namess[9])[2],')*(',JA2[i],')')
    }
  }
    if(checklist[7]==0)
   {
       Lam0=function(t){0}
   }
   if(checklist[8]==0)
   {
       Lam1=function(t){0}
   }
   if(checklist[9]==0)
   {
       Lam2=function(t){0}
   }

   odekernel=paste0('c(',dims[1],'\n,',dims[2],'\n,',dims[3],'\n,',dims[4],'\n,',
   dims2[1],'\n,',dims2[2],'\n,',dims2[3],'\n,',dims2[4],'\n,',
   body(namess[7])[2],'\n,',paste('(',body(namess[8])[2],')*x[5]'),'\n,',paste('(',body(namess[9])[2],')*x[6]'),')')
   }

   if(Jdist=='Normal')
   {
     prem =
     'mm1 = mu
     mm2 = (mu^2)+ (sig^2)
     mm3 = (mu^3)+ 3* (mu^1)*(sig^2)
     mm4 = (mu^4)+ 6* (mu^2)*(sig^2)+3*(sig^4)
     mm5 = (mu^5)+ 10*(mu^3)*(sig^2)+15*(mu^1)*(sig^4)
     mm6 = (mu^6)+ 15*(mu^4)*(sig^2)+45*(mu^2)*(sig^4)+15*(sig^6)
     mm7 = (mu^7)+ 21*(mu^5)*(sig^2)+105*(mu^3)*(sig^4)+105*mu*(sig^6)
     mm8 = (mu^8)+ 28*(mu^6)*(sig^2)+210*(mu^4)*(sig^4)+420*(mu^2)*(sig^6)+105*(sig^8)
     '
     prem=paste0('mu = ',body('Jmu')[2],'\n','sig = ', body('Jsig')[2],'\n',prem)

   }
      if(Jdist=='Exponential')
   {
     prem =
     'mm1=1*(mua^1)
	    mm2=2*(mua^2)
    	mm3=6*(mua^3)
    	mm4=24*(mua^4)
    	mm5=120*(mua^5)
    	mm6=720*(mua^6)
    	mm7=5040*(mua^7)
    	mm8=40320*(mua^8)
      '
     prem=paste0('mua = ',body('Jlam')[2],'\n',prem)

   }

    if(Jdist=='Gamma')
   {

     prem =
     'mm1=(alph*bet)
	    mm2=(alph+1)*alph*bet^2
    	mm3=(alph+2)*(alph+1)*alph*bet^3
    	mm4=(alph+3)*(alph+2)*(alph+1)*alph*bet^4
    	mm5=(alph+4)*(alph+3)*(alph+2)*(alph+1)*alph*bet^5
    	mm6=(alph+5)*(alph+4)*(alph+3)*(alph+2)*(alph+1)*alph*bet^6
    	mm7=(alph+6)*(alph+5)*(alph+4)*(alph+3)*(alph+2)*(alph+1)*alph*bet^7
    	mm8=(alph+7)*(alph+6)*(alph+5)*(alph+4)*(alph+3)*(alph+2)*(alph+1)*alph*bet^8
      '
     prem=paste0('alph = ',body('Jalpha')[2],'\n',prem)
     prem=paste0('bet  = ',body('Jbeta')[2],'\n',prem)
   }

   if(Jdist=='Laplace')
   {
     print('yeah')
     prem =
     'mm1=0.5*(+2*a^1*b^0)
      mm2=0.5*(+2*a^2*b^0+4*a^0*b^2)
     	mm3=0.5*(+2*a^3*b^0+12*a^1*b^2)
     	mm4=0.5*(+2*a^4*b^0+24*a^2*b^2+48*a^0*b^4)
     	mm5=0.5*(+2*a^5*b^0+40*a^3*b^2+240*a^1*b^4)
     	mm6=0.5*(+2*a^6*b^0+60*a^4*b^2+720*a^2*b^4+1440*a^0*b^6)
     	mm7=0.5*(+2*a^7*b^0+84*a^5*b^2+1680*a^3*b^4+10080*a^1*b^6)
     	mm8=0.5*(+2*a^8*b^0+112*a^6*b^2+3360*a^4*b^4+40320*a^2*b^6+80640*a^0*b^8)
      '
      prem=paste0('a = ',body('Ja')[2],'\n',prem)
      prem=paste0('b  = ',body('Jb')[2],'\n',prem)
      
   }
   if(Jdist=='Uniform')
   {
     prem =
     '
     mm1=(b+a)/2
     mm2=(a^2+a*b +b^2)/3
     mm3=1/4*(a+b)*(a^2+b^2)
     mm4=1/5*(a^4+a^3*b+a^2*b^2+a*b^3+b^4)
     mm5=1/6*(a^5+a^4*b+a^3*b^2+a^2*b^3+a^1*b^4+b^5)
     mm6=1/7*(a^6+a^5*b^1+a^4*b^2+a^3*b^3+a^2*b^4+a^1*b^5+b^6)
     mm7=1/8*(a^7+a^6*b^1+a^5*b^2+a^4*b^3+a^3*b^4+a^2*b^5+a^1*b^6+b^7)
     mm8=1/9*(a^8+a^7*b^1+a^6*b^2+a^5*b^3+a^4*b^4+a^3*b^5+a^2*b^6+a^1*b^7+b^8)
     '

      prem=paste0('a = ',body('Ja')[2],'\n',prem)
      prem=paste0('b  = ',body('Jb')[2],'\n',prem)
     
   }

   odekernel=paste0('{',prem,odekernel,'}')
   write(odekernel,'res.txt')
   ff <- function(x,t){}

   body(ff) = (parse(text=odekernel))
     t.alpha=
     c(0.000000000000000000000000000000000000000000000000000000000000,
       0.100000000000000000000000000000000000000000000000000000000000,
       0.539357840802981787532485197881302436857273449701009015505500,
       0.809036761204472681298727796821953655285910174551513523258250,
       0.309036761204472681298727796821953655285910174551513523258250,
       0.981074190219795268254879548310562080489056746118724882027805,
       0.833333333333333333333333333333333333333333333333333333333333,
       0.354017365856802376329264185948796742115824053807373968324184,
       0.882527661964732346425501486979669075182867844268052119663791,
       0.642615758240322548157075497020439535959501736363212695909875,
       0.357384241759677451842924502979560464040498263636787304090125,
       0.117472338035267653574498513020330924817132155731947880336209,
       0.833333333333333333333333333333333333333333333333333333333333,
       0.309036761204472681298727796821953655285910174551513523258250,
       0.539357840802981787532485197881302436857273449701009015505500,
       0.100000000000000000000000000000000000000000000000000000000000,
       1.00000000000000000000000000000000000000000000000000000000000)


       N=round((t-s)/delt)
       ttt=seq(s,t,by=delt)

       solver=function()
       {
         MM=matrix(0,2*TR.order+3,N)
         MA= c(rep(Xs^(1:TR.order),2),0,0,0)
         tt=s
         MM[,1]=MA
         ii=0
         while(ii < N)
        {
           ii=ii+1
           x0=MA
           fx0=ff(x0,tt)
           x1=x0+delt*(0.1*fx0)
           fx1=ff(x1,tt+t.alpha[2]*delt)
           x2=x0+delt*(-0.915176561375291*fx0  +1.45453440217827*fx1)
           fx2=ff(x2,tt+t.alpha[3]*delt)
           x3=x0+delt*( 0.202259190301118*fx0  +0.606777570903354*fx2)
           fx3=ff(x3,tt+t.alpha[4]*delt)
           x4=x0+delt*( 0.184024714708644*fx0  +0.197966831227192*fx2-0.0729547847313633*fx3)
           fx4=ff(x4,tt+t.alpha[5]*delt)
           x5=x0+delt*( 0.0879007340206681*fx0 +0.410459702520261*fx3+0.482713753678866*fx4)
           fx5=ff(x5,tt+t.alpha[6]*delt)
           x6=x0+delt*(0.085970050490246*fx0   +0.330885963040722*fx3+0.48966295730945*fx4-0.0731856375070851*fx5)
           fx6=ff(x6,tt+t.alpha[7]*delt)
           x7=x0+delt*(0.120930449125334*fx0   +0.260124675758296*fx4+0.0325402621549091*fx5-0.0595780211817361*fx6)
           fx7=ff(x7,tt+t.alpha[8]*delt)
           x8=x0+delt*(0.110854379580391*fx0   -0.0605761488255006*fx5+0.321763705601778*fx6+0.510485725608063*fx7)
           fx8=ff(x8,tt+t.alpha[9]*delt)
           x9=x0+delt*(0.112054414752879*fx0   -0.144942775902866*fx5-0.333269719096257*fx6+0.49926922955688*fx7+0.509504608929686*fx8)
           fx9=ff(x9,tt+t.alpha[10]*delt)
           x10=x0+delt*(0.113976783964186*fx0  -0.0768813364203357*fx5+0.239527360324391*fx6+0.397774662368095*fx7+0.0107558956873607*fx8-0.327769124164019*fx9)
           fx10=ff(x10,tt+t.alpha[11]*delt)
           x11=x0+delt*(0.0798314528280196*fx0 -0.0520329686800603*fx5-0.0576954146168549*fx6+0.194781915712104*fx7+0.145384923188325*fx8-0.0782942710351671*fx9-0.114503299361099*fx10)
           fx11=ff(x11,tt+t.alpha[12]*delt)
           x12=x0+delt*(0.985115610164857*fx0  +0.330885963040722*fx3+0.48966295730945*fx4-1.37896486574844*fx5-0.861164195027636*fx6+5.78428813637537*fx7+3.28807761985104*fx8-2.38633905093136*fx9-3.25479342483644*fx10-2.16343541686423*fx11)
           fx12=ff(x12,tt+t.alpha[13]*delt)
           x13=x0+delt*(0.895080295771633*fx0  +0.197966831227192*fx2-0.0729547847313633*fx3-0.851236239662008*fx5+0.398320112318533*fx6+3.63937263181036*fx7+1.5482287703983*fx8-2.12221714704054*fx9-1.58350398545326*fx10-1.71561608285936*fx11-0.0244036405750127*fx12)
           fx13=ff(x13,tt+t.alpha[14]*delt)
           x14=x0+delt*(-0.915176561375291*fx0+1.45453440217827*fx1+0*fx2+0*fx3-0.777333643644968*fx4+0*fx5-0.0910895662155176*fx6+0.0910895662155176*fx12+0.777333643644968*fx13)
           fx14=ff(x14,tt+t.alpha[15]*delt)
           x15=x0+delt*(0.1*fx0-0.157178665799771*fx2+0.157178665799771*fx14)
           fx15=ff(x15,tt+t.alpha[16]*delt)
           x16=x0+delt*(0.181781300700095*fx0+0.675*fx1+0.34275815984719*fx2+0*fx3+0.259111214548323*fx4-0.358278966717952*fx5-1.04594895940883*fx6+0.930327845415627*fx7+1.77950959431708*fx8+0.1*fx9-0.282547569539044*fx10-0.159327350119973*fx11-0.145515894647002*fx12-0.259111214548323*fx13-0.34275815984719*fx14-0.675*fx15)
           fx16=ff(x16,tt+t.alpha[17]*delt)

           MA=MA+(0.0333333333333333333333333333333333333333333333333333333333333*fx0
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
                +0.0333333333333333333333333333333333333333333333333333333333333*fx16)*delt

                tt=tt+delt
                MM[,ii]=MA
                setTxtProgressBar(pb,ii," "," ")
           }
           return(MM)
    }

  pb <- txtProgressBar(1,N,1,style = 1,width = 65)
  MM=solver()
  close(pb)

  if(eval.density)
  {
   DD=matrix(0,length(Xt),N)
      if(Dtype=='Gamma.A')
  {
    if(factorize)
    {
    if(factor.type=='Diffusion')
    {
     if(DTR.order==4)
    {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

       saddlep=function(xx,m)
       {
                 theta = (m[2]-m[1]^2)/m[1]
          kappa = m[1]/theta
          nu = dgamma(xx,shape=kappa,scale=theta)
          return(nu)
       }

       saddlep=function(xx,m)
       {
          theta = (m[2]-m[1]^2)/m[1]
          kappa = m[1]/theta
          nu = dgamma(xx,shape=kappa,scale=theta)
          return(nu)
       }
       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


      for(i in 1:N)
      {
         if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
         {
           DD[,i]=dnorm(Xt,MMM[1,i],sqrt(MMM[2,i]-MMM[1,i]^2))*exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
         }else
         {
           DD[,i]=saddlep(Xt,MM[1:4,i])
         }
         INDDDD[i] = indexs(MM[1:4,i])
      }

      DD=list(density=DD,MSH=0)

     }
    }
    if(factor.type=='Hawke')
    {
     if(DTR.order==4)
    {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^2)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^3)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^4)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

       saddlep=function(xx,m)
       {
          theta = (m[2]-m[1]^2)/m[1]
          kappa = m[1]/theta
          nu = dgamma(xx,shape=kappa,scale=theta)
          return(nu)
       }
       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


      for(i in 1:N)
      {
         if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
         {
           DD[,i]=exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])*Xt*(abs(Xt-MM[1+TR.order,i])<0.99*diff(Xt)[1])+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
         }else
         {
           DD[,i]=saddlep(Xt,MM[1:4,i])
         }
         INDDDD[i] = indexs(MM[1:4,i])
      }

      DD=list(density=DD,MSH=0)

     }
    }


    }


    if(!factorize)
    {

      if(DTR.order==4)
      {

      MMM = MM[(TR.order+1):(2*TR.order),]
      VVV = MM*0

      VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))


        saddlep=function(xx,m)
        {
          theta = (m[2]-m[1]^2)/m[1]
          kappa = m[1]/theta
          nu = dgamma(xx,shape=kappa,scale=theta)
          return(nu)
        }
        for(i in 1:N)
        {

            DD[,i]=saddlep(Xt,MM[1:4,i])

        }

        DD=list(density=DD,MSH=0)
      }
    }
    if(DTR.order==6)
    {
      stop('Incorrect input: 6th order saddlepoint not supported at present!')
    }
    if(DTR.order==8)
    {
      stop('Incorrect input: 8th order saddlepoint not supported at present!')
    }
  }

   if(Dtype=='Normal.A')
  {
    if(factorize)
    {
    if(factor.type=='Diffusion')
    {
     if(DTR.order==4)
    {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

       saddlep=function(xx,m)
       {
        k =m*0
        k[1]=                                          m[1]
        k[2]=                               m[2]-1*k[1]*m[1]
        k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
        k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]

        p=1/3*(3*(k[4]/6)*k[2] - ((k[3]/2)^2))/((k[4]/6)^2)
        q=1/27*(27*((k[4]/6)^2)*(k[1]-xx) - 9*(k[4]/6)*(k[3]/2)*k[2] + 2*((k[3]/2)^3))/((k[4]/6)^3)
        chk=(q^2)/4 + (p^3)/27
        th=-(k[3]/2)/(3*(k[4]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)

        K=k[1]*th+(k[2]*th^2)/2+(k[3]*th^3)/6 +(k[4]*th^4)/24
        K1=k[1]+(k[2]*th)+(k[3]*th^2)/2+(k[4]*th^3)/6
        K2=k[2]+(k[3]*th)+(k[4]*th^2)/2
        K3=k[3]+(k[4]*th)
        K4=k[4]
        return(1/sqrt(2*pi*(K2))*exp(K-th*K1))
       }

       saddlep=function(xx,m)
       {
        k =m*0
        k[1]=                                          m[1]
        k[2]=                               m[2]-1*k[1]*m[1]
        k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
        k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]

        p=1/3*(3*(k[4]/6)*k[2] - ((k[3]/2)^2))/((k[4]/6)^2)
        q=1/27*(27*((k[4]/6)^2)*(k[1]-xx) - 9*(k[4]/6)*(k[3]/2)*k[2] + 2*((k[3]/2)^3))/((k[4]/6)^3)
        chk=(q^2)/4 + (p^3)/27
        th=-(k[3]/2)/(3*(k[4]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)

        K=k[1]*th+(k[2]*th^2)/2+(k[3]*th^3)/6 +(k[4]*th^4)/24
        K1=k[1]+(k[2]*th)+(k[3]*th^2)/2+(k[4]*th^3)/6
        K2=k[2]+(k[3]*th)+(k[4]*th^2)/2
        K3=k[3]+(k[4]*th)
        K4=k[4]
        return(1/sqrt(2*pi*(K2))*exp(K-th*K1))
       }
       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


      for(i in 1:N)
      {
         if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
         {
           DD[,i]=dnorm(Xt,MMM[1,i],sqrt(MMM[2,i]-MMM[1,i]^2))*exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
         }else
         {
           DD[,i]=saddlep(Xt,MM[1:4,i])
         }
         INDDDD[i] = indexs(MM[1:4,i])
      }

      DD=list(density=DD,MSH=0)

     }
    }
        if(factor.type=='Hawke')
    {
     if(DTR.order==4)
    {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^2)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^3)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^4)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

       saddlep=function(xx,m)
       {
        k =m*0
        k[1]=                                          m[1]
        k[2]=                               m[2]-1*k[1]*m[1]
        k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
        k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]

        p=1/3*(3*(k[4]/6)*k[2] - ((k[3]/2)^2))/((k[4]/6)^2)
        q=1/27*(27*((k[4]/6)^2)*(k[1]-xx) - 9*(k[4]/6)*(k[3]/2)*k[2] + 2*((k[3]/2)^3))/((k[4]/6)^3)
        chk=(q^2)/4 + (p^3)/27
        th=-(k[3]/2)/(3*(k[4]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)

        K=k[1]*th+(k[2]*th^2)/2+(k[3]*th^3)/6 +(k[4]*th^4)/24
        K1=k[1]+(k[2]*th)+(k[3]*th^2)/2+(k[4]*th^3)/6
        K2=k[2]+(k[3]*th)+(k[4]*th^2)/2
        K3=k[3]+(k[4]*th)
        K4=k[4]
        return(1/sqrt(2*pi*(K2))*exp(K-th*K1))
       }
       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


      for(i in 1:N)
      {
         if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
         {
           DD[,i]=exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])*Xt*(abs(Xt-MM[1+TR.order,i])<0.99*diff(Xt)[1])+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
         }else
         {
           DD[,i]=saddlep(Xt,MM[1:4,i])
         }
         INDDDD[i] = indexs(MM[1:4,i])
      }

      DD=list(density=DD,MSH=0)

     }
    }


    }


    if(!factorize)
    {

      if(DTR.order==4)
      {

      MMM = MM[(TR.order+1):(2*TR.order),]
      VVV = MM*0

      VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))


        saddlep=function(xx,m)
        {
          k =m*0
          k[1]=                                          m[1]
          k[2]=                               m[2]-1*k[1]*m[1]

          return(dnorm(xx,k[1],sqrt(k[2])))
        }
        for(i in 1:N)
        {

            DD[,i]=saddlep(Xt,MM[1:4,i])

        }

        DD=list(density=DD,MSH=0)
      }
    }
    if(DTR.order==6)
    {
      stop('Incorrect input: 6th order saddlepoint not supported at present!')
    }
    if(DTR.order==8)
    {
      stop('Incorrect input: 8th order saddlepoint not supported at present!')
    }
  }
  if(Dtype=='Saddlepoint')
  {
    if(factorize)
    {
    if(factor.type=='Diffusion')
    {
     if(DTR.order==4)
    {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

       saddlep=function(xx,m)
       {
        k =m*0
        k[1]=                                          m[1]
        k[2]=                               m[2]-1*k[1]*m[1]
        k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
        k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]

        p=1/3*(3*(k[4]/6)*k[2] - ((k[3]/2)^2))/((k[4]/6)^2)
        q=1/27*(27*((k[4]/6)^2)*(k[1]-xx) - 9*(k[4]/6)*(k[3]/2)*k[2] + 2*((k[3]/2)^3))/((k[4]/6)^3)
        chk=(q^2)/4 + (p^3)/27
        th=-(k[3]/2)/(3*(k[4]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)

        K=k[1]*th+(k[2]*th^2)/2+(k[3]*th^3)/6 +(k[4]*th^4)/24
        K1=k[1]+(k[2]*th)+(k[3]*th^2)/2+(k[4]*th^3)/6
        K2=k[2]+(k[3]*th)+(k[4]*th^2)/2
        K3=k[3]+(k[4]*th)
        K4=k[4]
        return(1/sqrt(2*pi*(K2))*exp(K-th*K1))
       }
       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


      for(i in 1:N)
      {
         if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
         {
           DD[,i]=saddlep(Xt,MMM[1:4,i])*exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
         }else
         {
           DD[,i]=saddlep(Xt,MM[1:4,i])
         }
         INDDDD[i] = indexs(MM[1:4,i])
      }

      DD=list(density=DD,MSH=0)

     }
    }
        if(factor.type=='Hawke')
    {
     if(DTR.order==4)
    {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^2)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^3)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,]^4)/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

       saddlep=function(xx,m)
       {
        k =m*0
        k[1]=                                          m[1]
        k[2]=                               m[2]-1*k[1]*m[1]
        k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
        k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]

        p=1/3*(3*(k[4]/6)*k[2] - ((k[3]/2)^2))/((k[4]/6)^2)
        q=1/27*(27*((k[4]/6)^2)*(k[1]-xx) - 9*(k[4]/6)*(k[3]/2)*k[2] + 2*((k[3]/2)^3))/((k[4]/6)^3)
        chk=(q^2)/4 + (p^3)/27
        th=-(k[3]/2)/(3*(k[4]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)

        K=k[1]*th+(k[2]*th^2)/2+(k[3]*th^3)/6 +(k[4]*th^4)/24
        K1=k[1]+(k[2]*th)+(k[3]*th^2)/2+(k[4]*th^3)/6
        K2=k[2]+(k[3]*th)+(k[4]*th^2)/2
        K3=k[3]+(k[4]*th)
        K4=k[4]
        return(1/sqrt(2*pi*(K2))*exp(K-th*K1))
       }
       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


      for(i in 1:N)
      {
         if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
         {
           DD[,i]=exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])*Xt*(abs(Xt-MM[1+TR.order,i])<0.99*diff(Xt)[1])+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
         }else
         {
           DD[,i]=saddlep(Xt,MM[1:4,i])
         }
         INDDDD[i] = indexs(MM[1:4,i])
      }

      DD=list(density=DD,MSH=0)

     }
    }


    }


    if(!factorize)
    {

      if(DTR.order==4)
      {

      MMM = MM[(TR.order+1):(2*TR.order),]
      VVV = MM*0

      VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))


        saddlep=function(xx,m)
        {
          k =m*0
          k[1]=                                          m[1]
          k[2]=                               m[2]-1*k[1]*m[1]
          k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
          k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]

          p=1/3*(3*(k[4]/6)*k[2] - ((k[3]/2)^2))/((k[4]/6)^2)
          q=1/27*(27*((k[4]/6)^2)*(k[1]-xx) - 9*(k[4]/6)*(k[3]/2)*k[2] + 2*((k[3]/2)^3))/((k[4]/6)^3)
          chk=(q^2)/4 + (p^3)/27
          th=-(k[3]/2)/(3*(k[4]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)

          K=k[1]*th+(k[2]*th^2)/2+(k[3]*th^3)/6 +(k[4]*th^4)/24
          K1=k[1]+(k[2]*th)+(k[3]*th^2)/2+(k[4]*th^3)/6
          K2=k[2]+(k[3]*th)+(k[4]*th^2)/2
          K3=k[3]+(k[4]*th)
          K4=k[4]
          return(1/sqrt(2*pi*(K2))*exp(K-th*K1))
        }
        for(i in 1:N)
        {

            DD[,i]=saddlep(Xt,MM[1:4,i])

        }

        DD=list(density=DD,MSH=0)
      }
    }
    if(DTR.order==6)
    {
      stop('Incorrect input: 6th order saddlepoint not supported at present!')
    }
    if(DTR.order==8)
    {
      stop('Incorrect input: 8th order saddlepoint not supported at present!')
    }
  }
  if(Dtype=='Edgeworth')
  {
     edgeworth=function(x,m)
    {

        k =m*0
        k[1]=                                          m[1]
        k[2]=                               m[2]-1*k[1]*m[1]
        k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
        k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
        v3=k[3]/(sqrt(k[2])^3)
        v4=k[4]/(sqrt(k[2])^4)
        z=(x-k[1])/sqrt(k[2])
        h2=function(z){z^2-1}
        h3=function(z){z^3-3*z}
        h4=function(z){z^4-6*z^2+3}
        h5=function(z){z^5-10*z^3+15*z}
        return(pmax(dnorm(x,k[1],sqrt(k[2]))*(1+v3/6*h3(z)+v4/24*h4(z)),0))
    }

    if(factorize)
    {

    if(factor.type=='Diffusion')
    {
      if(DTR.order==4)
     {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))


       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


       for(i in 1:N)
       {
          if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
          {
            DD[,i]=edgeworth(Xt,MMM[1:4,i])*exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])+edgeworth(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
          }else
          {
            DD[,i]=edgeworth(Xt,MM[1:4,i])
          }
          INDDDD[i] = indexs(MM[1:4,i])
       }

       DD=list(density=DD,MSH=0)

      }
     }

     if(factor.type=='Hawke')
    {
      if(DTR.order==4)
     {

     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0
     p=exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])
     VVV[1,] = MM[1,]#-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-p*MM[1+TR.order,]^2)/(1-p)
     VVV[3,] = (MM[3,]-p*MM[1+TR.order,]^3)/(1-p)
     VVV[4,] = (MM[4,]-p*MM[1+TR.order,]^4)/(1-p)


       INDDDD = rep(0,N)
       indexs=function(m)
       {
         k =m*0
         k[1]=                                          m[1]
         k[2]=                               m[2]-1*k[1]*m[1]
         k[3]=                m[3]-1*k[1]*m[2]-2*k[2]*m[1]
         k[4]= m[4]-1*k[1]*m[3]-3*k[2]*m[2]-3*k[3]*m[1]
         skew = k[3]/(k[2]^(1.5))
         kurt = k[4]/(k[2]^2)
         return((skew^2+1)/kurt)
       }


       for(i in 1:N)
       {
          if(exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])>0.01)
          {
            DD[,i]=exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i])*Xt*(abs(Xt-MM[1,i])<0.005)+edgeworth(Xt,VVV[1:4,i])*(1-exp(-MM[2*TR.order+1,i]-MM[2*TR.order+2,i]-MM[2*TR.order+3,i]))#+dnorm(Xt,VVV[1,i],sqrt(VVV[2,i]-VVV[1,i]^2))*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))#+saddlep(Xt,VVV[1:4,i])*(1-exp(-MM[17,i]-MM[18,i]-MM[19,i]))
          }else
          {
            DD[,i]=edgeworth(Xt,MM[1:4,i])
          }
          INDDDD[i] = indexs(MM[1:4,i])
       }

       DD=list(density=DD,MSH=0)

      }
     }
    }
    if(!factorize)
    {
      if(DTR.order==4)
      {

      MMM = MM[(TR.order+1):(2*TR.order),]
      VVV = MM*0

      VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
      VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

        for(i in 1:N)
        {
            DD[,i]=edgeworth(Xt,MM[1:4,i])
        }

        DD=list(density=DD,MSH=0)
      }
    }
    if(DTR.order==6)
    {
      stop('Incorrect input: 6th order Edgeworth not supported at present!')
    }
    if(DTR.order==8)
    {
      stop('Incorrect input: 8th order Edgeworth not supported at present!')
    }
  }

   if(DTR.order==4)
  {
    pearson=function(x,k,type='Normal')
    {
      dx=diff(x)[1]
      u = k*0
      u[1,]=k[1,]
      u[2,]=k[2,]
      u[3,]=k[3,]
      u[4,]=k[4,]

      det=1*u[2,]*u[4,]+u[1,]*u[3,]*u[2,]+u[2,]*u[1,]*u[3,]-1*u[3,]*u[3,]-u[2,]*u[2,]*u[2,]-u[1,]*u[1,]*u[4,]

      b11=(u[2,]*u[4,]-u[3,]*u[3,])/det
      b12=(u[2,]*u[3,]-u[1,]*u[4,])/det
      b13=(u[1,]*u[3,]-u[2,]*u[2,])/det

      b21=(u[3,]*u[2,]-u[1,]*u[4,])/det
      b22=(1*u[4,]-u[2,]*u[2,])/det
      b23=(u[2,]*u[1,]-1*u[3,])/det

      b31=(u[1,]*u[3,]-u[2,]*u[2,])/det
      b32=(u[1,]*u[2,]-1*u[3,])/det
      b33=(1*u[2,]-u[1,]*u[1,])/det

      #Ntype
      if(type=='Normal')
      {
        V1=0
        V2=1
        V3=2*u[1,]
        betas1 =b11*V1+b12*V2+b13*V3
        betas2 =b21*V1+b22*V2+b23*V3
        betas3 =b31*V1+b32*V2+b33*V3

        thetas1=-betas1/1
        thetas2=-betas2/2
        thetas3=-betas3/3

        DD=matrix(0,length(x),length(u[1,]))
        MSH =matrix(0,P,length(u[1,]))
        MSH[1,] = lower
        MSH[P,] = upper
        for(i in 1:length(u[1,]))
        {
          y=seq(-1/(2*(u[1,i]-lower))*(sqrt(exp(2*alpha)+4*(u[1,i]-lower)^2)-exp(alpha)),-1/(2*(u[1,i]-upper))*(sqrt(exp(2*alpha)+4*(u[1,i]-upper)^2)-exp(alpha)),length=P)

          y=y[-c(1,length(y))]
          rho   = exp(alpha)*(1+y^2)/(1-(y)^2)^2
          tau    = exp(alpha)*(y)/(1-(y)^2)+u[1,i]
          MSH[2:(P-1),i] = tau
          dtemp  = sum(exp(thetas1[i]*tau+thetas2[i]*tau^2+thetas3[i]*tau^3)*rho*diff(y)[1])
          DD[,i] = exp(thetas1[i]*x+thetas2[i]*x^2+thetas3[i]*x^3)/dtemp
        }
        return(list(density=DD,MSH=MSH))

      }
      # Gtype
      if(type=='Gamma')
      {
        V1=1
        V2=2*u[1,]
        V3=3*u[2,]

        betas1 =b11*V1+b12*V2+b13*V3
        betas2 =b21*V1+b22*V2+b23*V3
        betas3 =b31*V1+b32*V2+b33*V3

        alphas  =1-betas1
        thetas1=-betas2/1
        thetas2=-betas3/2


        DD=matrix(0,length(x),length(u[1,]))
        MSH =matrix(0,P,length(u[1,]))
        MSH[1,] = lower
        MSH[P,] = upper
        for(i in 1:length(u[1,]))
        {
          y=seq(-1/(2*(u[1,i]-lower))*(sqrt(exp(2*alpha)+4*(u[1,i]-lower)^2)-exp(alpha)),-1/(2*(u[1,i]-upper))*(sqrt(exp(2*alpha)+4*(u[1,i]-upper)^2)-exp(alpha)),length=P)

          y=y[-c(1,length(y))]
          rho   = exp(alpha)*(1+y^2)/(1-(y)^2)^2
          tau    = exp(alpha)*(y)/(1-(y)^2)+u[1,i]
          MSH[2:(P-1),i] = tau
          dtemp  = sum(tau^(alphas[i]-1)*exp(thetas1[i]*tau+thetas2[i]*tau^2)*rho*diff(y)[1])
          DD[,i] = (x^(alphas[i]-1)*exp(thetas1[i]*x+thetas2[i]*x^2))/dtemp
        }
        return(list(density=DD,MSH=MSH))

      }

      # IGtype
      if(type=='InvGamma')
      {
        V1=2*u[1,]
        V2=3*u[2,]
        V3=4*u[3,]

        betas1 =b11*V1+b12*V2+b13*V3
        betas2 =b21*V1+b22*V2+b23*V3
        betas3 =b31*V1+b32*V2+b33*V3

        alphas  =-betas2
        nu     =betas1
        thetas =-betas3

        DD=matrix(0,length(x),length(u[1,]))
        MSH =matrix(0,P,length(u[1,]))
        MSH[1,] = lower
        MSH[P,] = upper
        for(i in 1:length(u[1,]))
        {
          y=seq(-1/(2*(u[1,i]-lower))*(sqrt(exp(2*alpha)+4*(u[1,i]-lower)^2)-exp(alpha)),-1/(2*(u[1,i]-upper))*(sqrt(exp(2*alpha)+4*(u[1,i]-upper)^2)-exp(alpha)),length=P)

          y=y[-c(1,length(y))]
          rho   = exp(alpha)*(1+y^2)/(1-(y)^2)^2
          tau    = exp(alpha)*(y)/(1-(y)^2)+u[1,i]
          MSH[2:(P-1),i] = tau
          dtemp  = sum((tau^alphas[i])*exp(nu[i]/tau+thetas[i]*tau)*rho*diff(y)[1])
          DD[,i] = (x^alphas[i])*exp(nu[i]/x+thetas[i]*x)/dtemp
        }
        return(list(density=DD,MSH=MSH))
      }

      # Btype
      if(type=='Beta')
      {
        V1=1-2*u[1,]
        V2=2*u[1,]-3*u[2,]
        V3=3*u[2,]-4*u[3,]

        betas1 =b11*V1+b12*V2+b13*V3
        betas2 =b21*V1+b22*V2+b23*V3
        betas3 =b31*V1+b32*V2+b33*V3

        alpha  =1-betas1
        nu     =1+(betas1+betas2+betas3)
        theta1 =betas3

        DD=matrix(0,length(x),length(u[1,]))
        for(i in 1:length(x))
        {
          DD[i,]=(x[i]^(alpha-1))*((1-x[i])^(nu-1))*exp(theta1*x[i])
        }
        for(i in 1:dim(u)[2])
        {
          DD[,i]=DD[,i]/sum(DD[,i]*dx)
        }
        return(list(density=DD,MSH=0))
      }

    }
  }
    if(DTR.order==8)
  {
    pearson=function(x,k,type='Normal')
    {
      dx=diff(x)[1]
      u1= k[1,]
      u2=k[2,]
      u3=k[3,]
      u4=k[4,]
      u5=k[5,]
      u6=k[6,]
      u7=k[7,]
      u8=k[8,]

      det=-u1*(u1*(u4*(u6*u8-u7*u7)-u5*(u5*u8-u6*u7)+u6*(u5*u7-u6*u6))-u3*(u2*(u6*u8-u7*u7)-u5*(u3*u8-u4*u7)+u6*(u3*u7-u4*u6))+u4*
                 (u2*(u5*u8-u6*u7)-u4*(u3*u8-u4*u7)+u6*(u3*u6-u4*u5))-u5*(u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)+u5*(u3*u6-u4*u5)))+u2*(u1*
                                                                                                                                     (u3*(u6*u8-u7*u7)-u5*(u4*u8-u5*u7)+u6*(u4*u7-u5*u6))-u2*(u2*(u6*u8-u7*u7)-u5*(u3*u8-u4*u7)+u6*(u3*u7-u4*u6))+u4*
                                                                                                                                     (u2*(u4*u8-u5*u7)-u3*(u3*u8-u4*u7)+(u3*u5-u4*u4)*u6)-u5*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u4*u6)+u5*(u3*u5-u4*u4)))-u3*(u1*
                                                                                                                                                                                                                                                         (u3*(u5*u8-u6*u7)-u4*(u4*u8-u5*u7)+u6*(u4*u6-u5*u5))-u2*(u2*(u5*u8-u6*u7)-u4*(u3*u8-u4*u7)+u6*(u3*u6-u4*u5))+u3*
                                                                                                                                                                                                                                                         (u2*(u4*u8-u5*u7)-u3*(u3*u8-u4*u7)+(u3*u5-u4*u4)*u6)-u5*(u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)+u4*(u3*u5-u4*u4)))+u2*
        (u4*(u6*u8-u7*u7)-u5*(u5*u8-u6*u7)+u6*(u5*u7-u6*u6))-u3*(u3*(u6*u8-u7*u7)-u5*(u4*u8-u5*u7)+u6*(u4*u7-u5*u6))+u4*
        (u3*(u5*u8-u6*u7)-u4*(u4*u8-u5*u7)+u6*(u4*u6-u5*u5))+u4*(u1*(u3*(u5*u7-u6*u6)-u4*(u4*u7-u5*u6)+u5*(u4*u6-u5*u5))-u2*
                                                                   (u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)+u5*(u3*u6-u4*u5))+u3*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u4*u6)+u5*(u3*u5-u4*u4))-u4*
                                                                   (u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)+u4*(u3*u5-u4*u4)))-u5*(u3*(u5*u7-u6*u6)-u4*(u4*u7-u5*u6)+u5*(u4*u6-u5*u5))

      b11= (u2*(u4*(u6*u8-u7*u7)-u5*(u5*u8-u6*u7)+u6*(u5*u7-u6*u6))-u3*(u3*(u6*u8-u7*u7)-u5*(u4*u8-u5*u7)+u6*(u4*u7-u5*u6))+u4*(u3*(u5*u8-u6*u7)-u4*(u4*u8-u5*u7)+u6*(u4*u6-u5*u5))-u5*(u3*(u5*u7-u6*u6)-u4*(u4*u7-u5*u6)+u5*(u4*u6-u5*u5)))/det
      b12= (-u1*(u4*(u6*u8-u7*u7)-u5*(u5*u8-u6*u7)+u6*(u5*u7-u6*u6))+u2*(u3*(u6*u8-u7*u7)-u5*(u4*u8-u5*u7)+u6*(u4*u7-u5*u6))-u3*(u3*(u5*u8-u6*u7)-u4*(u4*u8-u5*u7)+u6*(u4*u6-u5*u5))+u4*(u3*(u5*u7-u6*u6)-u4*(u4*u7-u5*u6)+u5*(u4*u6-u5*u5)))/det
      b13= (u1*(u3*(u6*u8-u7*u7)-u4*(u5*u8-u6*u7)+u5*(u5*u7-u6*u6))-u2*(u2*(u6*u8-u7*u7)-u4*(u4*u8-u5*u7)+u5*(u4*u7-u5*u6))+u3*(u2*(u5*u8-u6*u7)-u3*(u4*u8-u5*u7)+u5*(u4*u6-u5*u5))-u4*(u2*(u5*u7-u6*u6)-u3*(u4*u7-u5*u6)+u4*(u4*u6-u5*u5)))/det
      b14= (-u1*(u3*(u5*u8-u6*u7)-u4*(u4*u8-u6*u6)+u5*(u4*u7-u5*u6))+u2*(u2*(u5*u8-u6*u7)-u4*(u3*u8-u5*u6)+u5*(u3*u7-u5*u5))-u3*(u2*(u4*u8-u6*u6)-u3*(u3*u8-u5*u6)+u5*(u3*u6-u4*u5))+u4*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u5*u5)+u4*(u3*u6-u4*u5)))/det
      b15= (u1*(u3*(u5*u7-u6*u6)-u4*(u4*u7-u5*u6)+u5*(u4*u6-u5*u5))-u2*(u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)+u5*(u3*u6-u4*u5))+u3*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u4*u6)+u5*(u3*u5-u4*u4))-u4*(u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)+u4*(u3*u5-u4*u4)))/det

      b21= (-u1*(u4*(u6*u8-u7*u7)-u5*(u5*u8-u6*u7)+u6*(u5*u7-u6*u6))+u3*(u2*(u6*u8-u7*u7)-u5*(u3*u8-u4*u7)+u6*(u3*u7-u4*u6))-u4*(u2*(u5*u8-u6*u7)-u4*(u3*u8-u4*u7)+u6*(u3*u6-u4*u5))+u5*(u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)+u5*(u3*u6-u4*u5)))/det
      b22= (-u2*(u2*(u6*u8-u7*u7)-u5*(u3*u8-u4*u7)+u6*(u3*u7-u4*u6))+u3*(u2*(u5*u8-u6*u7)-u4*(u3*u8-u4*u7)+u6*(u3*u6-u4*u5))+u4*(u6*u8-u7*u7)-u5*(u5*u8-u6*u7)-u4*(u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)+u5*(u3*u6-u4*u5))+u6*(u5*u7-u6*u6))/det
      b23= (u2*(u1*(u6*u8-u7*u7)-u4*(u3*u8-u4*u7)+u5*(u3*u7-u4*u6))-u3*(u1*(u5*u8-u6*u7)-u3*(u3*u8-u4*u7)+u5*(u3*u6-u4*u5))-u3*(u6*u8-u7*u7)+u4*(u5*u8-u6*u7)+u4*(u1*(u5*u7-u6*u6)-u3*(u3*u7-u4*u6)+u4*(u3*u6-u4*u5))-u5*(u5*u7-u6*u6))/det
      b24= (-u2*(u1*(u5*u8-u6*u7)-u4*(u2*u8-u4*u6)+u5*(u2*u7-u4*u5))+u3*(u1*(u4*u8-u6*u6)-u3*(u2*u8-u4*u6)+u5*(u2*u6-u4*u4))+u3*(u5*u8-u6*u7)-u4*(u4*u8-u6*u6)-u4*(u1*(u4*u7-u5*u6)-u3*(u2*u7-u4*u5)+u4*(u2*u6-u4*u4))+u5*(u4*u7-u5*u6))/det
      b25= (u2*(u1*(u5*u7-u6*u6)-u4*(u2*u7-u3*u6)+u5*(u2*u6-u3*u5))-u3*(u1*(u4*u7-u5*u6)-u3*(u2*u7-u3*u6)+u5*(u2*u5-u3*u4))-u3*(u5*u7-u6*u6)+u4*(u4*u7-u5*u6)+u4*(u1*(u4*u6-u5*u5)-u3*(u2*u6-u3*u5)+u4*(u2*u5-u3*u4))-u5*(u4*u6-u5*u5))/det

      b31= (u1*(u3*(u6*u8-u7*u7)-u5*(u4*u8-u5*u7)+u6*(u4*u7-u5*u6))-u2*(u2*(u6*u8-u7*u7)-u5*(u3*u8-u4*u7)+u6*(u3*u7-u4*u6))+u4*(u2*(u4*u8-u5*u7)-u3*(u3*u8-u4*u7)+(u3*u5-u4*u4)*u6)-u5*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u4*u6)+u5*(u3*u5-u4*u4)))/det
      b32= (u1*(u2*(u6*u8-u7*u7)-u5*(u3*u8-u4*u7)+u6*(u3*u7-u4*u6))-u3*(u2*(u4*u8-u5*u7)-u3*(u3*u8-u4*u7)+(u3*u5-u4*u4)*u6)-u3*(u6*u8-u7*u7)+u5*(u4*u8-u5*u7)+u4*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u4*u6)+u5*(u3*u5-u4*u4))-u6*(u4*u7-u5*u6))/det
      b33= (-u1*(u1*(u6*u8-u7*u7)-u4*(u3*u8-u4*u7)+u5*(u3*u7-u4*u6))+u3*(u1*(u4*u8-u5*u7)-u2*(u3*u8-u4*u7)+u5*(u3*u5-u4*u4))+u2*(u6*u8-u7*u7)-u4*(u4*u8-u5*u7)-u4*(u1*(u4*u7-u5*u6)-u2*(u3*u7-u4*u6)+u4*(u3*u5-u4*u4))+u5*(u4*u7-u5*u6))/det
      b34= (u1*(u1*(u5*u8-u6*u7)-u4*(u2*u8-u4*u6)+u5*(u2*u7-u4*u5))-u3*(u1*(u3*u8-u5*u6)-u2*(u2*u8-u4*u6)+u5*(u2*u5-u3*u4))-u2*(u5*u8-u6*u7)+u4*(u3*u8-u5*u6)+u4*(u1*(u3*u7-u5*u5)-u2*(u2*u7-u4*u5)+u4*(u2*u5-u3*u4))-u5*(u3*u7-u5*u5))/det
      b35= (-u1*(u1*(u5*u7-u6*u6)-u4*(u2*u7-u3*u6)+u5*(u2*u6-u3*u5))+u3*(u1*(u3*u7-u4*u6)-u2*(u2*u7-u3*u6)+(u2*u4-u3*u3)*u5)+u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)-u4*(u1*(u3*u6-u4*u5)-u2*(u2*u6-u3*u5)+u4*(u2*u4-u3*u3))+u5*(u3*u6-u4*u5))/det

      b41= (-u1*(u3*(u5*u8-u6*u7)-u4*(u4*u8-u5*u7)+u6*(u4*u6-u5*u5))+u2*(u2*(u5*u8-u6*u7)-u4*(u3*u8-u4*u7)+u6*(u3*u6-u4*u5))-u3*(u2*(u4*u8-u5*u7)-u3*(u3*u8-u4*u7)+(u3*u5-u4*u4)*u6)+u5*(u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)+u4*(u3*u5-u4*u4)))/det
      b42= (-u1*(u2*(u5*u8-u6*u7)-u4*(u3*u8-u4*u7)+u6*(u3*u6-u4*u5))+u2*(u2*(u4*u8-u5*u7)-u3*(u3*u8-u4*u7)+(u3*u5-u4*u4)*u6)+u3*(u5*u8-u6*u7)-u4*(u4*u8-u5*u7)-u4*(u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)+u4*(u3*u5-u4*u4))+u6*(u4*u6-u5*u5))/det
      b43= (u1*(u1*(u5*u8-u6*u7)-u3*(u3*u8-u4*u7)+u5*(u3*u6-u4*u5))-u2*(u1*(u4*u8-u5*u7)-u2*(u3*u8-u4*u7)+u5*(u3*u5-u4*u4))-u2*(u5*u8-u6*u7)+u3*(u4*u8-u5*u7)+u4*(u1*(u4*u6-u5*u5)-u2*(u3*u6-u4*u5)+u3*(u3*u5-u4*u4))-u5*(u4*u6-u5*u5))/det
      b44= (-u1*(u1*(u4*u8-u6*u6)-u3*(u2*u8-u4*u6)+u5*(u2*u6-u4*u4))+u2*(u1*(u3*u8-u5*u6)-u2*(u2*u8-u4*u6)+u5*(u2*u5-u3*u4))+u2*(u4*u8-u6*u6)-u3*(u3*u8-u5*u6)-u4*(u1*(u3*u6-u4*u5)-u2*(u2*u6-u4*u4)+u3*(u2*u5-u3*u4))+u5*(u3*u6-u4*u5))/det
      b45= (u1*(u1*(u4*u7-u5*u6)-u3*(u2*u7-u3*u6)+u5*(u2*u5-u3*u4))-u2*(u1*(u3*u7-u4*u6)-u2*(u2*u7-u3*u6)+(u2*u4-u3*u3)*u5)-u2*(u4*u7-u5*u6)+u3*(u3*u7-u4*u6)+u4*(u1*(u3*u5-u4*u4)-u2*(u2*u5-u3*u4)+u3*(u2*u4-u3*u3))-u5*(u3*u5-u4*u4))/det

      b51= (u1*(u3*(u5*u7-u6*u6)-u4*(u4*u7-u5*u6)+u5*(u4*u6-u5*u5))-u2*(u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)+u5*(u3*u6-u4*u5))+u3*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u4*u6)+u5*(u3*u5-u4*u4))-u4*(u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)+u4*(u3*u5-u4*u4)))/det
      b52= (u1*(u2*(u5*u7-u6*u6)-u4*(u3*u7-u4*u6)+u5*(u3*u6-u4*u5))-u2*(u2*(u4*u7-u5*u6)-u3*(u3*u7-u4*u6)+u5*(u3*u5-u4*u4))-u3*(u5*u7-u6*u6)+u4*(u4*u7-u5*u6)+u3*(u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)+u4*(u3*u5-u4*u4))-u5*(u4*u6-u5*u5))/det
      b53= (-u1*(u1*(u5*u7-u6*u6)-u3*(u3*u7-u4*u6)+u4*(u3*u6-u4*u5))+u2*(u1*(u4*u7-u5*u6)-u2*(u3*u7-u4*u6)+u4*(u3*u5-u4*u4))+u2*(u5*u7-u6*u6)-u3*(u4*u7-u5*u6)-u3*(u1*(u4*u6-u5*u5)-u2*(u3*u6-u4*u5)+u3*(u3*u5-u4*u4))+u4*(u4*u6-u5*u5))/det
      b54= (u1*(u1*(u4*u7-u5*u6)-u3*(u2*u7-u4*u5)+u4*(u2*u6-u4*u4))-u2*(u1*(u3*u7-u5*u5)-u2*(u2*u7-u4*u5)+u4*(u2*u5-u3*u4))-u2*(u4*u7-u5*u6)+u3*(u3*u7-u5*u5)+u3*(u1*(u3*u6-u4*u5)-u2*(u2*u6-u4*u4)+u3*(u2*u5-u3*u4))-u4*(u3*u6-u4*u5))/det
      b55= (-u1*(u1*(u4*u6-u5*u5)-u3*(u2*u6-u3*u5)+u4*(u2*u5-u3*u4))+u2*(u1*(u3*u6-u4*u5)-u2*(u2*u6-u3*u5)+u4*(u2*u4-u3*u3))+u2*(u4*u6-u5*u5)-u3*(u3*u6-u4*u5)-u3*(u1*(u3*u5-u4*u4)-u2*(u2*u5-u3*u4)+u3*(u2*u4-u3*u3))+u4*(u3*u5-u4*u4))/det


      #Ntype
      if(type=='Normal')
      {
        V1=0
        V2=1
        V3=2*u1
        V4=3*u2
        V5=4*u3

        betas1 =b11*V1+b12*V2+b13*V3+b14*V4 +b15*V5
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4 +b25*V5
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4 +b35*V5
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4 +b45*V5
        betas5 =b51*V1+b52*V2+b53*V3+b54*V4 +b55*V5


        thetas1=-betas1/1
        thetas2=-betas2/2
        thetas3=-betas3/3
        thetas4=-betas4/4
        thetas5=-betas5/5

        DD=matrix(0,length(x),length(u1))
        MSH =matrix(0,P,length(u1))
        MSH[1,] = lower
        MSH[P,] = upper
        for(i in 1:length(u1))
        {
          y=seq(-1/(2*(u1[i]-lower))*(sqrt(exp(2*alpha)+4*(u1[i]-lower)^2)-exp(alpha)),-1/(2*(u1[i]-upper))*(sqrt(exp(2*alpha)+4*(u1[i]-upper)^2)-exp(alpha)),length=P)
          y=y[-c(1,length(y))]
          rho   = exp(alpha)*(1+y^2)/(1-(y)^2)^2
          tau    = exp(alpha)*(y)/(1-(y)^2)+u1[i]
          MSH[2:(P-1),i] = tau
          dtemp  = sum(exp(thetas1[i]*tau+thetas2[i]*tau^2+thetas3[i]*tau^3+thetas4[i]*tau^4+thetas5[i]*tau^5)*rho*diff(y)[1])
          DD[,i] = (exp(thetas1[i]*x+thetas2[i]*x^2+thetas3[i]*x^3+thetas4[i]*x^4+thetas5[i]*x^5))/dtemp
        }

        return(list(density=DD,MSH=MSH))


      }
      # Gtype
      if(type=='Gamma')
      {
        V1=1
        V2=2*u1
        V3=3*u2
        V4=4*u3
        V5=5*u4

        betas1 =b11*V1+b12*V2+b13*V3+b14*V4 +b15*V5
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4 +b25*V5
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4 +b35*V5
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4 +b45*V5
        betas5 =b51*V1+b52*V2+b53*V3+b54*V4 +b55*V5

        alphas  =1-betas1
        thetas1=-betas2/1
        thetas2=-betas3/2
        thetas3=-betas4/3
        thetas4=-betas5/4

        DD=matrix(0,length(x),length(u1))
        MSH =matrix(0,P,length(u1))
        MSH[1,] = lower
        MSH[P,] = upper
        for(i in 1:length(u1))
        {
          y=seq(-1/(2*(u1[i]-lower))*(sqrt(exp(2*alpha)+4*(u1[i]-lower)^2)-exp(alpha)),-1/(2*(u1[i]-upper))*(sqrt(exp(2*alpha)+4*(u1[i]-upper)^2)-exp(alpha)),length=P)

          y=y[-c(1,length(y))]
          rho   = exp(alpha)*(1+y^2)/(1-(y)^2)^2
          tau    = exp(alpha)*(y)/(1-(y)^2)+u1[i]
          MSH[2:(P-1),i] = tau
          dtemp  = sum(tau^(alphas[i]-1)*exp(thetas1[i]*tau+thetas2[i]*tau^2+thetas3[i]*tau^3+thetas4[i]*tau^4)*rho*diff(y)[1])
          DD[,i] = (x^(alphas[i]-1)*exp(thetas1[i]*x+thetas2[i]*x^2+thetas3[i]*x^3+thetas4[i]*x^4))/dtemp
        }
        return(list(density=DD,MSH=MSH))


      }

      # IGtype
      if(type=='InvGamma')
      {
        V1=2*u1
        V2=3*u2
        V3=4*u3
        V4=5*u4
        V5=6*u5

        betas1 =b11*V1+b12*V2+b13*V3+b14*V4 +b15*V5
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4 +b25*V5
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4 +b35*V5
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4 +b45*V5
        betas5 =b51*V1+b52*V2+b53*V3+b54*V4 +b55*V5

        alphas  =-betas2
        nu     =betas1
        theta1 =-betas3
        theta2 =-betas4/2
        theta3 =-betas5/3


        DD=matrix(0,length(x),length(u1))
        MSH =matrix(0,P,length(u1))
        MSH[1,] = lower
        MSH[P,] = upper
        for(i in 1:length(u1))
        {
          y=seq(-1/(2*(u1[i]-lower))*(sqrt(exp(2*alpha)+4*(u1[i]-lower)^2)-exp(alpha)),-1/(2*(u1[i]-upper))*(sqrt(exp(2*alpha)+4*(u1[i]-upper)^2)-exp(alpha)),length=P)
          y=y[-c(1,length(y))]
          rho   = exp(alpha)*(1+y^2)/(1-(y)^2)^2
          tau    = exp(alpha)*(y)/(1-(y)^2)+u1[i]
          MSH[2:(P-1),i] = tau
          dtemp  = sum((tau^alphas[i])*exp(nu[i]/tau+theta1[i]*tau+theta2[i]*tau^2+theta3[i]*tau^3)*rho*diff(y)[1])
          DD[,i] = ((x^alphas[i])*exp(nu[i]/x+theta1[i]*x+theta2[i]*x^2+theta3[i]*x^3))/dtemp
        }

        return(list(density=DD,MSH=MSH))
      }

      # Btype
      if(type=='Beta')
      {
        V1=1-2*u1
        V2=2*u1-3*u2
        V3=3*u2-4*u3
        V4=4*u3-5*u4
        V5=5*u4-6*u5

        betas1 =b11*V1+b12*V2+b13*V3+b14*V4 +b15*V5
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4 +b25*V5
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4 +b35*V5
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4 +b45*V5
        betas5 =b51*V1+b52*V2+b53*V3+b54*V4 +b55*V5

        alpha  =1-betas1
        nu     =1+(betas1+betas2+betas3+betas4+betas5)
        theta1 =betas3+betas4+betas5
        theta2 =(betas4+betas5)/2
        theta3 =betas5/3

        DD=matrix(0,length(x),length(u1))
        for(i in 1:length(x))
        {
          DD[i,]=(x[i]^(alpha-1))*((1-x[i])^(nu-1))*exp(theta1*x[i]+theta2*x[i]^2+theta3*x[i]^3)
        }
        for(i in 1:length(u1))
        {
          DD[,i]=DD[,i]/sum(DD[,i]*dx)
        }
        return(list(density=DD,MSH=0))
      }

    }
  }
  if(Dtype=='Normal')
  {
     MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[5,] = (MM[5,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[5+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[6,] = (MM[6,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[6+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[7,] = (MM[7,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[7+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[8,] = (MM[8,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[8+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

    if(factorize)
    {
      prb=exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])
      d1=pearson(Xt,MMM,'Normal')
      #print(d1)
      d2=pearson(Xt,VVV,'Normal')
      DD=d1$density
      #print(dim(DD))
      for(i in 1:dim(d1$density)[2])
      {
        DD[,i]=d1$density[,i]*prb[i]+d2$density[,i]*(1-prb[i])
      }
      DD=list(density=DD)
    }
    if(!factorize)
    {
      DD=pearson(Xt,MMM,'Normal')
    }
  }
  if(Dtype=='Gamma')
  {
    if(min(Xt)<=0){stop('Incorrect input: Xt must be strictly positive for Dtype = Gamma. ')}
          MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

      VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[5,] = (MM[5,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[5+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[6,] = (MM[6,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[6+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[7,] = (MM[7,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[7+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[8,] = (MM[8,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[8+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

    if(factorize)
    {
      prb=exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])
      d1=pearson(Xt,MMM,'Gamma')
      #print(d1)
      d2=pearson(Xt,VVV,'Gamma')
      DD=d1$density
      #print(dim(DD))
      for(i in 1:dim(d1$density)[2])
      {
        DD[,i]=d1$density[,i]*prb[i]+d2$density[,i]*(1-prb[i])
      }
      DD=list(density=DD)
    }
    if(!factorize)
    {
      DD=pearson(Xt,MMM,'Gamma')
    }
  }
  if(Dtype=='InvGamma')
  {

    if(min(Xt)<=0){stop('Incorrect input: Xt must be strictly positive for Dtype = InvGamma. ')}
    DD=pearson(Xt,MM,'InvGamma')
  }
  if(Dtype=='Beta')
  {
    if(min(Xt)<=0){stop('Incorrect input: Xt must in [0,1] Dtype = Beta. Try seq(delta,1-delta,delta) for delta>0.')}
         MMM = MM[(TR.order+1):(2*TR.order),]
     VVV = MM*0

     VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[5,] = (MM[5,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[5+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[6,] = (MM[6,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[6+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[7,] = (MM[7,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[7+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
     VVV[8,] = (MM[8,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[8+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))

    if(factorize)
    {
      prb=exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])
      d1=pearson(Xt,MMM,'Beta')
      #print(d1)
      d2=pearson(Xt,VVV,'Beta')
      DD=d1$density
      #print(dim(DD))
      for(i in 1:dim(d1$density)[2])
      {
        DD[,i]=d1$density[,i]*prb[i]+d2$density[,i]*(1-prb[i])
      }
      DD=list(density=DD)
    }
    if(!factorize)
    {
      DD=pearson(Xt,MMM,'Beta')

    }
  }
}
  }
   if(!eval.density)
  {
    DD=list(density=NULL,MSH=NULL);
    MMM = MM[(TR.order+1):(2*TR.order),]
    VVV = MM*0

    VVV[1,] = (MM[1,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[1+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
    VVV[2,] = (MM[2,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[2+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
    VVV[3,] = (MM[3,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[3+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
    VVV[4,] = (MM[4,]-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,])*MM[4+TR.order,])/(1-exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]))
  }

  # Fix this!!!
  if(TR.order==8)
  {
  u=MM[1:8+TR.order,]
  k=u*0
  k[1,]=                                                                                                      u[1,]
  k[2,]=                                                                                         u[2,]-1*k[1,]*u[1,]
  k[3,]=                                                                           u[3,]-1*k[1,]*u[2,]-2*k[2,]*u[1,]
  k[4,]=                                                             u[4,]-1*k[1,]*u[3,]-3*k[2,]*u[2,]-3*k[3,]*u[1,]
  k[5,]=                                               u[5,]-1*k[1,]*u[4,]-4*k[2,]*u[3,]-6*k[3,]*u[2,]-4*k[4,]*u[1,]
  k[6,]=                               u[6,]-1*k[1,]*u[5,]-5*k[2,]*u[4,]-10*k[3,]*u[3,]-10*k[4,]*u[2,]-5*k[5,]*u[1,]
  k[7,]=                u[7,]-1*k[1,]*u[6,]-6*k[2,]*u[5,]-15*k[3,]*u[4,]-20*k[4,]*u[3,]-15*k[5,]*u[2,]-6*k[6,]*u[1,]
  k[8,]= u[8,]-1*k[1,]*u[7,]-7*k[2,]*u[6,]-21*k[3,]*u[5,]-35*k[4,]*u[4,]-35*k[5,]*u[3,]-21*k[6,]*u[2,]-7*k[7,]*u[1,]
  }
  if(TR.order==4)
  {
    u=MM[1:4+TR.order,]
    k=u*0
    k[1,]=                                                                                                      u[1,]
    k[2,]=                                                                                         u[2,]-1*k[1,]*u[1,]
    k[3,]=                                                                           u[3,]-1*k[1,]*u[2,]-2*k[2,]*u[1,]
    k[4,]=                                                             u[4,]-1*k[1,]*u[3,]-3*k[2,]*u[2,]-3*k[3,]*u[1,]
  }


  ret= list(density=DD$density,Xt=Xt,time=seq(s,t,delt)[-1],cumulants=k,moments=MM,mesh=DD$MSH,zero_jump_prob =exp(-MM[2*TR.order+1,]-MM[2*TR.order+2,]-MM[2*TR.order+3,]),excess_moments =VVV)
  class(ret) = 'JGQD.density'
  return(ret)

}
