GQD.density<-function(Xs,Xt,s,t,delt=1/100,Dtype='Saddlepoint',Trunc=c(4,4),P=100,alpha=0,lower=0,upper=50,print.output=TRUE,eval.density=TRUE)
{
   lookin=1
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
  Dtypes =c('Saddlepoint','Normal','Gamma','InvGamma','Beta')
  Dindex = which(Dtypes==Dtype)
  
  b1 = '\n==============================================================================\n'
  b2 = '==============================================================================\n'
  warn=c(
     'Input 1: Argument {lower} must be of length 1.\n'
    ,'Input 2: Argument {upper} must be of length 1.\n'
    ,'Input 3: {upper} must be > {lower}.\n'
    ,'Input 4: Range [lower,upper] must be strictly positive for Dtype Gamma or InvGamma.\n'
    ,'Input 5: Dtype cannot be Beta for observations not in (0,1).\n'
    ,'Input 6: P must be >= 10.\n'
    ,'Input 7: Dtype has to be one of Saddlepoint, Normal, Gamma, InvGamma or Beta.\n'
    ,'Input 8: Trunc[2] must be <= Trunc[1].\n'
    ,'Input 9: Large {delt} may result in poor approximations.\n'
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
    ,'Input 20: Cumulant truncation (Trunc[1]) must be one of 4, 6 or 8.\n'
    ,'Input 21: Density truncation (Trunc[2]) must be one of 4, 6 or 8.\n'
    ,'Input 22: Starting time {s} cannot be greater than {t}.\n'
  )
  
   warntrue=rep(FALSE,40)
   
  if(length(Xt)<2){warntrue[14] =TRUE}#{stop(paste0(b1,warn[14],b2))}
  if(length(Xs)!=1){warntrue[15] =TRUE}#{stop(paste0(b1,warn[15],b2))}
  if(!is.vector(Xt)){warntrue[16] =TRUE}#{stop(paste0(b1,warn[16],b2))}
  if(sum(Dindex)==0){warntrue[7] =TRUE}#{stop(paste0(b1,warn[7],b2))}
  if(length(lower)>1){warntrue[1] =TRUE}#{stop(paste0(b1,warn[1],b2))}
  if(length(upper)>1){warntrue[2] =TRUE}#{stop(paste0(b1,warn[2],b2))}
  if(upper<=lower){warntrue[3] =TRUE}#{stop(paste0(b1,warn[3],b2))}
  if((Dindex==3)|(Dindex==4)){if(lower[1]<=0){warntrue[4] =TRUE}}
  if(Dindex==5){if(any(Xt<=0)|any(Xt>=1)){warntrue[5] =TRUE}}
  if(P<10){warntrue[6] =TRUE}#{stop(paste0(b1,warn[6],b2))}
  if(length(P)!=1){warntrue[11] =TRUE}#{stop(paste0(b1,warn[11],b2))}
  if(length(alpha)!=1){warntrue[18] =TRUE}#{stop(paste0(b1,warn[18],b2))}
  if(Trunc[2]>Trunc[1]){warntrue[8] =TRUE}#{stop(paste0(b1,warn[8],b2))}
  if(length(Trunc)!=2){warntrue[19] =TRUE}#{stop(paste0(b1,warn[19],b2))}
  if(sum(TR.order==c(4,6,8))!=1){warntrue[20] =TRUE}#{stop(paste0(b1,warn[20],b2))}
  if(sum(DTR.order==c(4,6,8))!=1){warntrue[21] =TRUE}#{stop(paste0(b1,warn[21],b2))}
  #if(delt>=1/5){stop(paste0(b1,warn[9],b2))}
  if(delt>=1){warntrue[10] =TRUE}#{stop(paste0(b1,warn[10],b2))}
  if(t<s){warntrue[22] =TRUE}#{stop(paste0(b1,warn[22],b2))}

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

  pow=function(x,p)
  {
    x^p
  }
  prod=function(a,b){a*b}
  
  if(sum(TR.order==c(4,6,8))!=1)
  {
    stop('Incorrect input: TR.order must be of order 4, 6 or 8!')
  }
  if(sum(Dtype==c('Saddlepoint','Normal','Gamma','InvGamma','Beta'))!=1)
  {
    stop('Incorrect input: Dtype must be one of Saddle, Normal, Gamma, InvGamma or Beta!')
  }
  function.list=objects(pos=lookin)
  coef.terms=c('G0','G1','G2','Q0','Q1','Q2')
  coef.terms2=paste0(coef.terms,'(t)')
  checklist=rep(0,6)
  for(i in 1:6)
  {
    checklist[i]= sum(function.list== coef.terms[i])
  }
  checklist
  coef.index=which(checklist==1)
  coef.index
  
  if(TR.order==4)
  {
    M4 = rbind(
      c('1','(y[1])','(1*y[2]+1*y[1]*y[1])'                              ,' ',' ',' '),
      c(' ','(2*y[2])','(2*y[3]+4*y[1]*y[2])'                            ,'1','(   y[1])','(y[2]+y[1]*y[1])'),
      c(' ','(3*y[3])','(3*y[4]+6*y[1]*y[3]+6*y[2]*y[2])'                ,' ','( 3*y[2])','(3*y[3]+6*y[1]*y[2])'),
      c(' ','(4*y[4])','(      +8*y[1]*y[4]+24*y[2]*y[3])'               ,' ','( 6*y[3])','(6*y[4]+12*y[1]*y[3]+12*y[2]*y[2])'))
    
    
    dims=rep('',4)
    for(i in coef.index)
    {
      for(j in 1:4)
      {
        if(M4[j,i]!=' ')
        {
          dims[j] = paste0(dims[j],'+(',body(coef.terms[i])[2],')*',M4[j,i])
        }
      }
    }
    if(dims[1]==''){dims[1]='0'}
    
    res=paste0('c(',dims[1],',',dims[2],',',dims[3],',',dims[4],')')
    if(sum(checklist[c(3,5,6)])==0)
    {
      res=paste0('c(',dims[1],',',dims[2],')')
      TR.order=2
      Dtype='Saddlepoint'
    }
    ff=function(y,t){}
    body(ff)= parse(text = res)
  }
  if(TR.order==6)
  {
    M6 = rbind(
      c('1','(y[1])','(1*y[2]+1*y[1]*y[1])'                              ,' ',' ',' '),
      c(' ','(2*y[2])','(2*y[3]+4*y[1]*y[2])'                            ,'1','(   y[1])','(y[2]+y[1]*y[1])'),
      c(' ','(3*y[3])','(3*y[4]+6*y[1]*y[3]+6*y[2]*y[2])'                ,' ','( 3*y[2])','(3*y[3]+6*y[1]*y[2])'),
      c(' ','(4*y[4])','(4*y[5]+8*y[1]*y[4]+24*y[2]*y[3])'               ,' ','( 6*y[3])','(6*y[4]+12*y[1]*y[3]+12*y[2]*y[2])'),
      c(' ','(5*y[5])','(5*y[6]+10*y[1]*y[5]+40*y[2]*y[4]+30*y[3]*y[3])' ,' ','(10*y[4])','(10*y[5]+20*y[1]*y[4]+60*y[2]*y[3])'),
      c(' ','(6*y[6])','(      +12*y[1]*y[6]+60*y[2]*y[5]+120*y[3]*y[4])',' ','(15*y[5])','(15*y[6]+30*y[1]*y[5]+120*y[2]*y[4]+90*y[3]*y[3])'))
    
    
    dims=rep('',6)
    for(i in coef.index)
    {
      for(j in 1:6)
      {
        if(M6[j,i]!=' ')
        {
          dims[j] = paste0(dims[j],'+(',body(coef.terms[i])[2],')*',M6[j,i])
        }
      }
    }
    if(dims[1]==''){dims[1]='0'}
    res=paste0('c(',dims[1],',',dims[2],',',dims[3],',',dims[4],',',dims[5],',',dims[6],')')
    if(sum(checklist[c(3,5,6)])==0)
    {
      res=paste0('c(',dims[1],',',dims[2],')')
      TR.order=2
      Dtype='Saddlepoint'
    }
    ff=function(y,t){}
    body(ff)= parse(text = res)
  }
  if(TR.order==8)
  {
    M8 = rbind(
      c('1','  (y[1])','(1*y[2]+1*y[1]*y[1])'                                               ,' ',' ',' '),
      c(' ','(2*y[2])','(2*y[3]+4*y[1]*y[2])'                                               ,'1','(   y[1])','(y[2]+y[1]*y[1])'),
      c(' ','(3*y[3])','(3*y[4]+6*y[1]*y[3]+6*y[2]*y[2])'                                   ,' ','( 3*y[2])','(3*y[3]+6*y[1]*y[2])'),
      c(' ','(4*y[4])','(4*y[5]+8*y[1]*y[4]+24*y[2]*y[3])'                                  ,' ','( 6*y[3])','(6*y[4]+12*y[1]*y[3]+12*y[2]*y[2])'),
      c(' ','(5*y[5])','(5*y[6]+10*y[1]*y[5]+40*y[2]*y[4]+30*y[3]*y[3])'                    ,' ','(10*y[4])','(10*y[5]+20*y[1]*y[4]+60*y[2]*y[3])'),
      c(' ','(6*y[6])','(6*y[7]+12*y[1]*y[6]+60*y[2]*y[5]+120*y[3]*y[4])'                   ,' ','(15*y[5])','(15*y[6]+30*y[1]*y[5]+120*y[2]*y[4]+90*y[3]*y[3])'),
      c(' ','(7*y[7])','(7*y[8]+2*7*y[1]*y[7]+2*42*y[2]*y[6]+2*105*y[3]*y[5]+140*y[4]*y[4])',' ','(21*y[6])','(21*y[7]+42*y[1]*y[6]+210*y[2]*y[5]+420*y[3]*y[4])'),
      c(' ','(8*y[8])','(      +16*y[1]*y[8]+112*y[2]*y[7]+2*168*y[3]*y[6]+2*280*y[5]*y[4])',' ','(28*y[7])','(28*y[8]+56*y[1]*y[7]+336*y[2]*y[6]+840*y[3]*y[5]+560*y[4]*y[4])'))
    
    dims=rep('',8)
    for(i in coef.index)
    {
      for(j in 1:8)
      {
        if(M8[j,i]!=' ')
        {
          dims[j] = paste0(dims[j],'+(',body(coef.terms[i])[2],')*',M8[j,i])
        }
      }
    }
    if(dims[1]==''){dims[1]='0'}
    res=paste0('c(',dims[1],',',dims[2],',',dims[3],',',dims[4],',',dims[5],',',dims[6],',',dims[7],',',dims[8],')')
    if(sum(checklist[c(3,5,6)])==0)
    {
      res=paste0('c(',dims[1],',',dims[2],')')
      TR.order=2
      Dtype='Saddlepoint'
    }
    ff=function(y,t){}
    body(ff)= parse(text = res)
  }
  
  Dsuggest=F
  if(Dsuggest)
  {
    Dtype='skip'
    if(sum(checklist[c(3,5,6)])==0)
    {
      Dtype='Saddlepoint'
      print('Dsuggest: Saddlepoint is exact. Not running addtional approximations.',quote=F)
    }
  }
  
  #==============================================================================
  #                           Interface Module
  #==============================================================================
  
  namess4=c('G0','G1','G2','Q0','Q1','Q2')
  trim <- function (x) gsub("([[:space:]])", "", x)
  
  for(i in 1:6)
  {
    if(sum(function.list==namess4[i]))
    {
      namess4[i]=paste0(namess4[i],' : ',trim(body(namess4[i])[2]))
    }
  }
  namess4=matrix(namess4,length(namess4),1)
  
  dinfo = c('Density approx. : ',
            'P               : ',
            'alpha           : ',
            'Trunc. Order    : ',
            'Dens.  Order    : ')
  dinfo[1] =paste0(dinfo[1],Dtype)
  if(Dindex!=1)
  {
    dinfo[2] =paste0(dinfo[2],P)
    dinfo[3] =paste0(dinfo[3],alpha)
  }
  dinfo[4] =paste0(dinfo[4],Trunc[1])
  dinfo[5] =paste0(dinfo[5],Trunc[2])
  
  buffer0=c('================================================================')
  buffer1=c('----------------------------------------------------------------')
  buffer2=c('................................................................')
  buffer3=c('...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ')
  buffer4=c('_____________________ Drift Coefficients _______________________')
  buffer5=c('___________________ Diffusion Coefficients _____________________')
  buffer6=c('__________________ Distribution Approximant ____________________')
  
  type.sol ="                 Generalized Quadratic Diffusion (GQD) "
  
  Info=c(buffer0,type.sol,buffer0,buffer4,namess4[1:3],buffer5,namess4[4:6],buffer6,dinfo)
  Info=data.frame(matrix(Info,length(Info),1))
  colnames(Info)=''
  if(print.output)
  {
  print(Info,row.names = FALSE,right=F)
  }
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
    MM=matrix(0,TR.order,N)
    MA= c(Xs,rep(0,TR.order-1))
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
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  if(DTR.order==4)
  {
    pearson=function(x,k,type='Normal')
    {
      dx=diff(x)[1]
      u = k*0
      u[1,]=                                          k[1,]
      u[2,]=                            k[2,]+1*k[1,]*u[1,]
      u[3,]=              k[3,]+1*k[1,]*u[2,]+2*k[2,]*u[1,]
      u[4,]=k[4,]+1*k[1,]*u[3,]+3*k[2,]*u[2,]+3*k[3,]*u[1,]
      
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
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  
  if(DTR.order==6)
  {
    pearson=function(x,k,type='Normal')
    {
      dx=diff(x)[1]
      u = k*0
      u[1,]=                                          k[1,]
      u[2,]=                            k[2,]+1*k[1,]*u[1,]
      u[3,]=              k[3,]+1*k[1,]*u[2,]+2*k[2,]*u[1,]
      u[4,]=k[4,]+1*k[1,]*u[3,]+3*k[2,]*u[2,]+3*k[3,]*u[1,]
      
      u[1,]=                                                                        k[1,]
      u[2,]=                                                          k[2,]+1*k[1,]*u[1,]
      u[3,]=                                            k[3,]+1*k[1,]*u[2,]+2*k[2,]*u[1,]
      u[4,]=                              k[4,]+1*k[1,]*u[3,]+3*k[2,]*u[2,]+3*k[3,]*u[1,]
      u[5,]=                k[5,]+1*k[1,]*u[4,]+4*k[2,]*u[3,]+6*k[3,]*u[2,]+4*k[4,]*u[1,]
      u[6,]=k[6,]+1*k[1,]*u[5,]+5*k[2,]*u[4,]+10*k[3,]*u[3,]+10*k[4,]*u[2,]+5*k[5,]*u[1,]
      
      det=-u[1,]*(u[1,]*(u[4,]*u[6,]-u[5,]*u[5,])-u[3,]*(u[2,]*u[6,]-u[3,]*u[5,])+u[4,]*(u[2,]*u[5,]-u[3,]*u[4,]))+u[2,]*(u[1,]*(u[3,]*u[6,]-u[4,]*u[5,])-u[2,]*(u[2,]*u[6,]-u[3,]*u[5,])+u[4,]*(u[2,]*u[4,]-u[3,]*u[3,]))+u[2,]*(u[4,]*u[6,]-u[5,]*u[5,])-u[3,]*(u[3,]*u[6,]-u[4,]*u[5,])-u[3,]*(u[1,]*(u[3,]*u[5,]-u[4,]*u[4,])-u[2,]*(u[2,]*u[5,]-u[3,]*u[4,])+u[3,]*(u[2,]*u[4,]-u[3,]*u[3,]))+u[4,]*(u[3,]*u[5,]-u[4,]*u[4,])
      
      b11=(u[2,]*(u[4,]*u[6,]-u[5,]*u[5,])-u[3,]*(u[3,]*u[6,]-u[4,]*u[5,])+u[4,]*(u[3,]*u[5,]-u[4,]*u[4,]))/det
      b12=(-u[1,]*(u[4,]*u[6,]-u[5,]*u[5,])+u[2,]*(u[3,]*u[6,]-u[4,]*u[5,])-u[3,]*(u[3,]*u[5,]-u[4,]*u[4,]))/det
      b13=(u[1,]*(u[3,]*u[6,]-u[4,]*u[5,])-u[2,]*(u[2,]*u[6,]-u[4,]*u[4,])+u[3,]*(u[2,]*u[5,]-u[3,]*u[4,]))/det
      b14=(-u[1,]*(u[3,]*u[5,]-u[4,]*u[4,])+u[2,]*(u[2,]*u[5,]-u[3,]*u[4,])-u[3,]*(u[2,]*u[4,]-u[3,]*u[3,]))/det
      
      b21=(-u[1,]*(u[4,]*u[6,]-u[5,]*u[5,])+u[3,]*(u[2,]*u[6,]-u[3,]*u[5,])-u[4,]*(u[2,]*u[5,]-u[3,]*u[4,]))/det
      b22=(-u[2,]*(u[2,]*u[6,]-u[3,]*u[5,])+u[4,]*u[6,]-u[5,]*u[5,]+u[3,]*(u[2,]*u[5,]-u[3,]*u[4,]))/det
      b23=(u[2,]*(u[1,]*u[6,]-u[3,]*u[4,])-u[3,]*u[6,]-u[3,]*(u[1,]*u[5,]-u[3,]*u[3,])+u[4,]*u[5,])/det
      b24=(-u[2,]*(u[1,]*u[5,]-u[2,]*u[4,])+u[3,]*u[5,]-u[4,]*u[4,]+u[3,]*(u[1,]*u[4,]-u[2,]*u[3,]))/det
      
      b31=(u[1,]*(u[3,]*u[6,]-u[4,]*u[5,])-u[2,]*(u[2,]*u[6,]-u[3,]*u[5,])+u[4,]*(u[2,]*u[4,]-u[3,]*u[3,]))/det
      b32=(u[1,]*(u[2,]*u[6,]-u[3,]*u[5,])-u[3,]*u[6,]+u[4,]*u[5,]-u[3,]*(u[2,]*u[4,]-u[3,]*u[3,]))/det
      b33=(-u[1,]*(u[1,]*u[6,]-u[3,]*u[4,])+u[2,]*u[6,]-u[4,]*u[4,]+u[3,]*(u[1,]*u[4,]-u[2,]*u[3,]))/det
      b34=(u[1,]*(u[1,]*u[5,]-u[2,]*u[4,])-u[2,]*u[5,]+u[3,]*u[4,]-u[3,]*(u[1,]*u[3,]-u[2,]*u[2,]))/det
      
      b41=(-u[1,]*(u[3,]*u[5,]-u[4,]*u[4,])+u[2,]*(u[2,]*u[5,]-u[3,]*u[4,])-u[3,]*(u[2,]*u[4,]-u[3,]*u[3,]))/det
      b42=(-u[1,]*(u[2,]*u[5,]-u[3,]*u[4,])+u[3,]*u[5,]-u[4,]*u[4,]+u[2,]*(u[2,]*u[4,]-u[3,]*u[3,]))/det
      b43=(u[1,]*(u[1,]*u[5,]-u[3,]*u[3,])-u[2,]*u[5,]-u[2,]*(u[1,]*u[4,]-u[2,]*u[3,])+u[3,]*u[4,])/det
      b44=(-u[1,]*(u[1,]*u[4,]-u[2,]*u[3,])+u[2,]*u[4,]-u[3,]*u[3,]+u[2,]*(u[1,]*u[3,]-u[2,]*u[2,]))/det
      
      #Ntype
      if(type=='Normal')
      {
        V1=0
        V2=1
        V3=2*u[1,]
        V4=3*u[2,]
        
        betas1 =b11*V1+b12*V2+b13*V3+b14*V4
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4
        
        thetas1=-betas1/1
        thetas2=-betas2/2
        thetas3=-betas3/3
        thetas4=-betas4/4
        
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
          dtemp  = sum(exp(thetas1[i]*tau+thetas2[i]*tau^2+thetas3[i]*tau^3+thetas4[i]*tau^4)*rho*diff(y)[1])
          DD[,i] = exp(thetas1[i]*x+thetas2[i]*x^2+thetas3[i]*x^3+thetas4[i]*x^4)/dtemp
        }
        return(list(density=DD,MSH=MSH))
      }
      # Gtype
      if(type=='Gamma')
      {
        V1=1
        V2=2*u[1,]
        V3=3*u[2,]
        V4=4*u[3,]
        
        betas1 =b11*V1+b12*V2+b13*V3+b14*V4
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4
        
        alphas  =1-betas1
        thetas1=-betas2/1
        thetas2=-betas3/2
        thetas3=-betas4/3
        
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
          dtemp  = sum(tau^(alphas[i]-1)*exp(thetas1[i]*tau+thetas2[i]*tau^2+thetas3[i]*tau^3)*rho*diff(y)[1])
          DD[,i] = (x^(alphas[i]-1)*exp(thetas1[i]*x+thetas2[i]*x^2+thetas3[i]*x^3))/dtemp
        }
        return(list(density=DD,MSH=MSH))
        
        
      }
      
      # IGtype
      if(type=='InvGamma')
      {
        V1=2*u[1,]
        V2=3*u[2,]
        V3=4*u[3,]
        V4=5*u[4,]
        
        betas1 =b11*V1+b12*V2+b13*V3+b14*V4
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4
        
        alphas  =-betas2
        nu     =betas1
        theta1 =-betas3
        theta2 =-betas4/2
        
        DD=matrix(0,length(x),length(u[1,]))
        MSH=matrix(0,P,length(u[1,]))
        MSH[1,] = lower
        MSH[P,] = upper
        for(i in 1:length(u[1,]))
        {
          y=seq(-1/(2*(u[1,i]-lower))*(sqrt(exp(2*alpha)+4*(u[1,i]-lower)^2)-exp(alpha)),-1/(2*(u[1,i]-upper))*(sqrt(exp(2*alpha)+4*(u[1,i]-upper)^2)-exp(alpha)),length=P)
          y=y[-c(1,length(y))]
          rho   = exp(alpha)*(1+y^2)/(1-(y)^2)^2
          tau    = exp(alpha)*(y)/(1-(y)^2)+u[1,i]
          MSH[2:(P-1),i] = tau
          dtemp  = sum((tau^alphas[i])*exp(nu[i]/tau+theta1[i]*tau+theta2[i]*tau^2)*rho*diff(y)[1])
          DD[,i] = (x^alphas[i])*exp(nu[i]/x+theta1[i]*x+theta2[i]*x^2)/dtemp
        }
        return(list(density=DD,MSH=MSH))
      }
      
      # Btype
      if(type=='Beta')
      {
        V1=1-2*u[1,]
        V2=2*u[1,]-3*u[2,]
        V3=3*u[2,]-4*u[3,]
        V4=4*u[3,]-5*u[4,]
        
        betas1 =b11*V1+b12*V2+b13*V3+b14*V4
        betas2 =b21*V1+b22*V2+b23*V3+b24*V4
        betas3 =b31*V1+b32*V2+b33*V3+b34*V4
        betas4 =b41*V1+b42*V2+b43*V3+b44*V4
        
        alpha  =1-betas1
        nu     =1+(betas1+betas2+betas3+betas4)
        theta1 =betas3+betas4
        theta2 =betas4/2
        
        DD=matrix(0,length(x),length(u[1,]))
        for(i in 1:length(x))
        {
          DD[i,]=(x[i]^(alpha-1))*((1-x[i])^(nu-1))*exp(theta1*x[i]+theta2*x[i]^2)
        }
        for(i in 1:dim(u)[2])
        {
          DD[,i]=DD[,i]/sum(DD[,i]*dx)
        }
        return(list(density=DD,MSH=0))
      }
      
    }
  }
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  
  if(DTR.order==8)
  {
    pearson=function(x,k,type='Normal')
    {
      dx=diff(x)[1]
      u1=                                                                                  k[1,]
      u2=                                                                       k[2,]+1*k[1,]*u1
      u3=                                                            k[3,]+1*k[1,]*u2+2*k[2,]*u1
      u4=                                                 k[4,]+1*k[1,]*u3+3*k[2,]*u2+3*k[3,]*u1
      u5=                                      k[5,]+1*k[1,]*u4+4*k[2,]*u3+6*k[3,]*u2+4*k[4,]*u1
      u6=                         k[6,]+1*k[1,]*u5+5*k[2,]*u4+10*k[3,]*u3+10*k[4,]*u2+5*k[5,]*u1
      u7=             k[7,]+1*k[1,]*u6+6*k[2,]*u5+15*k[3,]*u4+20*k[4,]*u3+15*k[5,]*u2+6*k[6,]*u1
      u8= k[8,]+1*k[1,]*u7+7*k[2,]*u6+21*k[3,]*u5+35*k[4,]*u4+35*k[5,]*u3+21*k[6,]*u2+7*k[7,]*u1
      
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
  DD=matrix(0,length(Xt),N)
  if(Dtype=='Saddlepoint')
  {
    if(sum(checklist[c(3,5,6)])==0)
    {
      for(i in 1:length(Xt))
      {
        DD[i,]=dnorm(Xt[i],MM[1,],sd=sqrt(MM[2,]))
      }
      ret=list(Xt=Xt,time=ttt[-1],density=DD,cumulants=MM)
      class(ret) = 'GQD.density'
      return(ret)
    }
    if(DTR.order==4)
    {
      MSH = 0
      for(i in 1:length(Xt))
      {
        p=1/3*(3*(MM[4,]/6)*MM[2,] - ((MM[3,]/2)^2))/((MM[4,]/6)^2)
        q=1/27*(27*((MM[4,]/6)^2)*(MM[1,]-Xt[i]) - 9*(MM[4,]/6)*(MM[3,]/2)*MM[2,] + 2*((MM[3,]/2)^3))/((MM[4,]/6)^3)
        chk=(q^2)/4 + (p^3)/27
        th=-(MM[3,]/2)/(3*(MM[4,]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)
        
        K=MM[1,]*th+(MM[2,]*th^2)/2+(MM[3,]*th^3)/6 +(MM[4,]*th^4)/24
        K1=MM[1,]+(MM[2,]*th)+(MM[3,]*th^2)/2+(MM[4,]*th^3)/6
        K2=MM[2,]+(MM[3,]*th)+(MM[4,]*th^2)/2
        K3=MM[3,]+(MM[4,]*th)
        K4=MM[4,]
        DD[i,]=1/sqrt(2*pi*(K2))*exp(K-th*K1)#*(1+K4/(K2^2)/8-5*(K3^2)/(K2^3)/24)
      }
      DD=list(density=DD,MSH=0)
    }
    if(DTR.order==6)
    {
      #for(i in 1:length(Xt))
      #{
      #   p=1/3*(3*(MM[4,]/6)*MM[2,] - ((MM[3,]/2)^2))/((MM[4,]/6)^2)
      #   q=1/27*(27*((MM[4,]/6)^2)*(MM[1,]-Xt[i]) - 9*(MM[4,]/6)*(MM[3,]/2)*MM[2,] + 2*((MM[3,]/2)^3))/((MM[4,]/6)^3)
      #   chk=(q^2)/4 + (p^3)/27
      #   th=-(MM[3,]/2)/(3*(MM[4,]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)
      #
      #   K=MM[1,]*th+(MM[2,]*th^2)/2+(MM[3,]*th^3)/6 +(MM[4,]*th^4)/24
      #   K1=MM[1,]+(MM[2,]*th)+(MM[3,]*th^2)/2+(MM[4,]*th^3)/6
      #   K2=MM[2,]+(MM[3,]*th)+(MM[4,]*th^2)/2
      #   K3=MM[3,]+(MM[4,]*th)
      #   K4=MM[4,]
      #   DD[i,]=1/sqrt(2*pi*(K2))*exp(K-th*K1)#*(1+K4/(K2^2)/8-5*(K3^2)/(K2^3)/24)
      #}
      stop('Incorrect input: 6th order saddlepoint not supported at present!')
    }
    if(DTR.order==6)
    {
      #for(i in 1:length(Xt))
      #{
      #   p=1/3*(3*(MM[4,]/6)*MM[2,] - ((MM[3,]/2)^2))/((MM[4,]/6)^2)
      #   q=1/27*(27*((MM[4,]/6)^2)*(MM[1,]-Xt[i]) - 9*(MM[4,]/6)*(MM[3,]/2)*MM[2,] + 2*((MM[3,]/2)^3))/((MM[4,]/6)^3)
      #   chk=(q^2)/4 + (p^3)/27
      #   th=-(MM[3,]/2)/(3*(MM[4,]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)
      #
      #   K=MM[1,]*th+(MM[2,]*th^2)/2+(MM[3,]*th^3)/6 +(MM[4,]*th^4)/24
      #   K1=MM[1,]+(MM[2,]*th)+(MM[3,]*th^2)/2+(MM[4,]*th^3)/6
      #   K2=MM[2,]+(MM[3,]*th)+(MM[4,]*th^2)/2
      #   K3=MM[3,]+(MM[4,]*th)
      #   K4=MM[4,]
      #   DD[i,]=1/sqrt(2*pi*(K2))*exp(K-th*K1)#*(1+K4/(K2^2)/8-5*(K3^2)/(K2^3)/24)
      #}
      stop('Incorrect input: 8th order saddlepoint not supported at present!')
    }
  }
  if(Dtype=='Normal')
  {
    DD=pearson(Xt,MM,'Normal')
  }
  if(Dtype=='Gamma')
  {
    if(min(Xt)<=0){stop('Incorrect input: Xt must be strictly positive for Dtype = Gamma. ')}
    DD=pearson(Xt,MM,'Gamma')
  }
  if(Dtype=='InvGamma')
  {
    
    if(min(Xt)<=0){stop('Incorrect input: Xt must be strictly positive for Dtype = InvGamma. ')}
    DD=pearson(Xt,MM,'InvGamma')
  }
  if(Dtype=='Beta')
  {
    if(min(Xt)<=0){stop('Incorrect input: Xt must in [0,1] Dtype = Beta. Try seq(delta,1-delta,delta) for delta>0.')}
    DD=pearson(Xt,MM,'Beta')
  }
  
  # if(Dsuggest)
  #{
  #
  # if(sum(checklist[c(3,5,6)])==0)
  # {
  #    for(i in 1:length(Xt))
  #    {
  #        DD[i,]=dnorm(Xt[i],MM[1,],sd=sqrt(MM[2,]))
  #    }
  # }
  # if(sum(checklist[c(3,5,6)])!=0)
  # {
  #  for(i in 1:length(Xt))
  #  {
  #   p=1/3*(3*(MM[4,]/6)*MM[2,] - ((MM[3,]/2)^2))/((MM[4,]/6)^2)
  #    q=1/27*(27*((MM[4,]/6)^2)*(MM[1,]-Xt[i]) - 9*(MM[4,]/6)*(MM[3,]/2)*MM[2,] + 2*((MM[3,]/2)^3))/((MM[4,]/6)^3)
  #    chk=(q^2)/4 + (p^3)/27
  #    th=-(MM[3,]/2)/(3*(MM[4,]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)
  #
  #    K=MM[1,]*th+(MM[2,]*th^2)/2+(MM[3,]*th^3)/6 +(MM[4,]*th^4)/24
  #    K1=MM[1,]+(MM[2,]*th)+(MM[3,]*th^2)/2+(MM[4,]*th^3)/6
  #    K2=MM[2,]+(MM[3,]*th)+(MM[4,]*th^2)/2
  #    K3=MM[3,]+(MM[4,]*th)
  #    K4=MM[4,]
  #    DD[i,]=1/sqrt(2*pi*(K2))*exp(K-th*K1)#*(1+K4/(K2^2)/8-5*(K3^2)/(K2^3)/24)
  #  }
  #  DD1 = DD
  #  DD2=pearson(Xt,MM,'Normal')
  #  ifelse(min(Xt)<=0,
  #  {
  #     warning('Dsuggest: min(Xt)<=0. Modified Xt for Dtype = Gamma.')
  #     Xtmod = seq(diff(Xt)[1],max(Xt),diff(Xt)[1])
  #     DD3=pearson(Xtmod,MM,'Gamma')
  #  },
  #  {
  #     Xtmod=Xt
  #     DD3=pearson(Xtmod,MM,'Gamma')
  #  })
  #   ifelse(min(Xt)<=0,
  #  {
  #     warning('Dsuggest: min(Xt)<=0. Modified Xt for Dtype = InvGamma.')
  #     Xtmod = seq(diff(Xt)[1],max(Xt),diff(Xt)[1])
  #     DD4=pearson(Xtmod,MM,'InvGamma')
  #  },
  #  {
  #     Xtmod = Xt
  #     DD4=pearson(Xtmod,MM,'InvGamma')
  #  })
  #   ifelse((min(Xt)<=0)|(min(Xt)>=1),
  #  {
  #     warning('Dsuggest: Xt not in [0,1]. Modified Xt for Dtype = Beta.')
  #     Xtmod2 = seq(diff(Xt)[1],1-diff(Xt)[1],diff(Xt)[1])
  #     DD5=pearson(Xtmod2,MM,'Beta')
  #  },
  #  {
  #     Xtmod2=Xt
  #     DD5=pearson(Xtmod2,MM,'Beta')
  #  })
  #
  #  DD =list(Saddle=DD1,Pearson.Normal=DD2,Pearson.Gamma=DD3,Pearson.InvGamma=DD4,Pearson.Beta = DD5)
  #  Xt =list(Xt,Xt,Xtmod,Xtmod,Xtmod2)
  #
  #  calc.moment=function(x,D)
  #  {
  #     mom1 = 0*D[1,]
  #     dx=diff(x)[1]
  #     for(i in 1:length(x))
  #      {
  #        mom1=mom1+x[i]*D[i,]*dx
  #     }
  #     return((MM[1,]-mom1)^2)
  #  }
  #  deviation  = matrix(0,5,dim(MM)[2])
  #  for(i in 1:5)
  #  {
  #      deviation[i,]=calc.moment(Xt[[i]],DD[[i]])
  #  }
  
  #par(mfrow=c(1,2))
  #plot(1,1,type='n',xlim=range(ttt[-1]),ylim=c(0,1))
  #for(i in 1:4)
  #{
  #   lines(deviation[i,]~ttt[-1],col=i+1)
  #}
  #  mins=apply(deviation,2,which.min)
  #  plot(mins~ttt[-1],type='s',axes=F,main='Min Deviation',ylab='Type',xlab='Time(t)',ylim=c(0.75,5.25))
  #  axis(2,at=1:5,labels = names(DD),cex.axis=0.8)
  #  axis(1)
  
  #  counts=c(sum(mins==1),sum(mins==2),sum(mins==3),sum(mins==4),sum(mins==5))
  #  best=which.max(counts)
  #  ind=rep('',5)
  #  ind[best]='[=]'
  #  res=cbind(ind)
  #  rownames(res)=names(DD)
  #  data.frame(res)
  #  colnames(res) ='Suggested Density'
  #  print(res,right=T,quote=F)
  # }
  #}
  }
  if(TR.order==4)
  {
    u=MM*0
    k=MM
    u[1,]=                                                                                                      k[1,] 
    u[2,]=                                                                                        k[2,]+1*k[1,]*u[1,] 
    u[3,]=                                                                          k[3,]+1*k[1,]*u[2,]+2*k[2,]*u[1,] 
    u[4,]=                                                            k[4,]+1*k[1,]*u[3,]+3*k[2,]*u[2,]+3*k[3,]*u[1,]  
  }
  if(TR.order==6)
  {
    u=MM*0
    k=MM
    u[1,]=                                                                                                      k[1,] 
    u[2,]=                                                                                        k[2,]+1*k[1,]*u[1,] 
    u[3,]=                                                                          k[3,]+1*k[1,]*u[2,]+2*k[2,]*u[1,] 
    u[4,]=                                                            k[4,]+1*k[1,]*u[3,]+3*k[2,]*u[2,]+3*k[3,]*u[1,] 
    u[5,]=                                              k[5,]+1*k[1,]*u[4,]+4*k[2,]*u[3,]+6*k[3,]*u[2,]+4*k[4,]*u[1,] 
    u[6,]=                              k[6,]+1*k[1,]*u[5,]+5*k[2,]*u[4,]+10*k[3,]*u[3,]+10*k[4,]*u[2,]+5*k[5,]*u[1,]
  }
  if(TR.order==8)
  {
    u=MM*0
    k=MM
    u[1,]=                                                                                                      k[1,] 
    u[2,]=                                                                                        k[2,]+1*k[1,]*u[1,] 
    u[3,]=                                                                          k[3,]+1*k[1,]*u[2,]+2*k[2,]*u[1,] 
    u[4,]=                                                            k[4,]+1*k[1,]*u[3,]+3*k[2,]*u[2,]+3*k[3,]*u[1,] 
    u[5,]=                                              k[5,]+1*k[1,]*u[4,]+4*k[2,]*u[3,]+6*k[3,]*u[2,]+4*k[4,]*u[1,] 
    u[6,]=                              k[6,]+1*k[1,]*u[5,]+5*k[2,]*u[4,]+10*k[3,]*u[3,]+10*k[4,]*u[2,]+5*k[5,]*u[1,] 
    u[7,]=               k[7,]+1*k[1,]*u[6,]+6*k[2,]*u[5,]+15*k[3,]*u[4,]+20*k[4,]*u[3,]+15*k[5,]*u[2,]+6*k[6,]*u[1,] 
    u[8,]=k[8,]+1*k[1,]*u[7,]+7*k[2,]*u[6,]+21*k[3,]*u[5,]+35*k[4,]*u[4,]+35*k[5,]*u[3,]+21*k[6,]*u[2,]+7*k[7,]*u[1,] 
  }

  if(!eval.density)
  {
      DD=list(density=NULL,MSH=NULL);
  }
  ret=list(density=DD$density,Xt=Xt,time=ttt[-1],cumulants=MM,moments=u,mesh=DD$MSH)
  class(ret) = 'GQD.density'
  return(ret)
  }