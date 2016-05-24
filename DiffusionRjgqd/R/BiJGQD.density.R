

 BiJGQD.density=function(Xs,Ys,Xt=seq(1,18,1/2),Yt=seq(1,18,1/2),s,t,delt,Dtype='Saddlepoint',Jdist='MVNormal',Jtype='Add',print.output=TRUE,eval.density=TRUE)
 {

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
      a00=function(t){10+sin(2*pi*t)}
      a10=function(t){-1}
      c10=function(t){1}

      b00=function(t){10+cos(2*pi*t)}
      b01=function(t){-1}
      f01=function(t){2}

      model=BiGQD.density(10,10,seq(5,15,1/5),seq(5,15,1/5),0,5,delt=1/100)
      GQD.plot(model)
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
      c00=function(t){2+sin(2*pi*t)}
      f00=function(t){2+sin(2*pi*t)}
      model=BiGQD.density(10,10,seq(5,15,1/5),seq(5,15,1/5),0,3,delt=1/100)
      GQD.plot(model)
      --------------------------------------------------------------------------------
      '
      check=T
    }
    return(list(check=check,txt=txt))
  }
  check_for=check_for_model()
  if(check_for[[1]]){stop(check_for[[2]])}

  # Warning Module
  Dtypes =c('Saddlepoint','Edgeworth')
  Dindex = which(Dtypes==Dtype)
  b1 = '\n==============================================================================\n'
  b2 = '==============================================================================\n'
  warn=c(
    '1. Input: Dtype has to be one of "Saddlepoint" or "Edgeworth".\n'
    ,'2. Input: Argument {delt} must be < 1.\n'
    ,'3. Input: length(delt)!=1.\n'
    ,'4. Input: length(Xt)=length(Yt) must be > 1.\n'
    ,'5. Input: Arguments {Xs,Ys} must be of length 1.\n'
    ,'6. Input: Arguments {Xt,Yt} must be of type vector!.\n'
    ,'7. Input: Starting time {s} cannot be greater than {t}.\n'
    ,'8. Input: length(Xt)!=length(Yt).\n'
    ,'9. Input: Arguments {s,t} must be of length 1.\n'
    ,'10. Input: {Jdist} has to be of type "MVNormal".\n'
    ,'11. Input: {Jtype} has to be of type "Add" or "Mult".\n'
  )

  JDtypes=c('Normal','MVNormal')
  JDindex = which(JDtypes==Jdist)
  warntrue=rep(FALSE,20)
  if(length(Xt)<2){warntrue[4] =TRUE}#{stop(paste0(b1,warn[4],b2))}
  if(length(Xs)!=1){warntrue[5] =TRUE}#{stop(paste0(b1,warn[5],b2))}
  if(length(Yt)<2){warntrue[4] =TRUE}#{stop(paste0(b1,warn[4],b2))}
  if(length(Ys)!=1){warntrue[5] =TRUE}#{stop(paste0(b1,warn[5],b2))}
  if(!is.vector(Xt)){warntrue[6] =TRUE}#{stop(paste0(b1,warn[6],b2))}
  if(!is.vector(Yt)){warntrue[6] =TRUE}#{stop(paste0(b1,warn[6],b2))}
  if(sum(Dindex)==0){warntrue[1] =TRUE}#{stop(paste0(b1,warn[1],b2))}
  if(delt>=1){warntrue[2] =TRUE}#{stop(paste0(b1,warn[2],b2))}
  if(length(delt)>1){warntrue[3] =TRUE}#{stop(paste0(b1,warn[3],b2))}
  if(t<s){warntrue[7] =TRUE}#{stop(paste0(b1,warn[7],b2))}
  if(length(Xt)!=length(Yt)){warntrue[8] =TRUE}#{stop(paste0(b1,warn[8],b2))}
  if(length(t)!=1){warntrue[9] =TRUE}#{stop(paste0(b1,warn[9],b2))}
  if(length(s)!=1){warntrue[9] =TRUE}#{stop(paste0(b1,warn[9],b2))}
  if(sum(JDindex)==0){warntrue[10] =TRUE}
  if((Jtype!='Add')&&(Jtype!='Mult')){warntrue[11] =TRUE}
  if(eval.density)
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

  nnn=length(Xt)
  X1=matrix(outer(Xt,rep(1,nnn)),1,nnn*nnn,byrow=T)
  X2=matrix(t(outer(Yt,rep(1,nnn))),1,nnn*nnn,byrow=T)

  pow=function(x,p)
  {
    x^p
  }
  prod=function(a,b){a*b}

   X1=matrix(outer(Xt,rep(1,nnn)),1,nnn*nnn,byrow=T)
   X2=matrix(t(outer(Yt,rep(1,nnn))),1,nnn*nnn,byrow=T)




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

  Mstar2 =
  c("+1*mv10*m10"
   ,"+2*mv10*m20+1*mv20*m10"
   ,"+3*mv10*m30+3*mv20*m20+1*mv30*m10"
   ,"+4*mv10*m40+6*mv20*m30+4*mv30*m20+1*mv40*m10"
   ,"+1*mv01*m10"
   ,"+2*mv01*m11+1*mv02*m10"
   ,"+3*mv01*m12+3*mv02*m11+1*mv03*m10"
   ,"+4*mv01*m13+6*mv02*m12+4*mv03*m11+1*mv04*m10"
   ,"mv10*m11+mv01*m20+mv11*m10"
   ,"mv10*m12+2*mv01*m21+2*mv11*m11+mv02*m20+mv12*m10"
   ,"2*mv10*m21+mv20*m11+mv01*m30+2*mv11*m20+mv21*m10"
   ,"2*mv10*m22+mv20*m12+2*mv01*m31+4*mv11*m21+2*mv21*m11+mv02*m30+2*mv12*m20+mv22*m10"
   ,"mv10*m13+3*mv01*m22+3*mv11*m12+3*mv02*m21+3*mv12*m11+mv03*m20+mv13*m10"
   ,"3*mv10*m31+3*mv20*m21+mv30*m11+mv01*m40+3*mv11*m30+3*mv21*m20+mv31*m10")

   Mstar3 =
  c("+1*mv10*m01"
   ,"+2*mv10*m11+1*mv20*m01"
   ,"+3*mv10*m21+3*mv20*m11+1*mv30*m01"
   ,"+4*mv10*m31+6*mv20*m21+4*mv30*m10+1*mv40*m01"
   ,"+1*mv01*m01"
   ,"+2*mv01*m02+1*mv02*m01"
   ,"+3*mv01*m03+3*mv02*m02+1*mv03*m01"
   ,"+4*mv01*m04+6*mv02*m03+4*mv03*m02+1*mv04*m01"
   ,"mv10*m02+mv01*m11+mv11*m01"
   ,"mv10*m03+2*mv01*m12+2*mv11*m02+mv02*m11+mv12*m01"
   ,"2*mv10*m12+mv20*m02+mv01*m21+2*mv11*m11+mv21*m01"
   ,"2*mv10*m13+mv20*m03+2*mv01*m22+4*mv11*m12+2*mv21*m02+mv02*m21+2*mv12*m11+mv22*m01"
   ,"mv10*m04+3*mv01*m13+3*mv11*m03+3*mv02*m12+3*mv12*m02+mv03*m11+mv13*m01"
   ,"3*mv10*m22+3*mv20*m12+mv30*m02+mv01*m31+3*mv11*m21+3*mv21*m11+mv31*m01")
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

     Mstar2 =
  c("m20*(mv10)"
   ,"m30*(mv20+2*mv10)"
   ,"m40*(mv30+3*mv20+3*mv10)"
   ,"m50*(mv40+4*mv30+6*mv20+4*mv10)"
   ,"m11*(mv01)"
   ,"m12*(mv02+2*mv01)"
   ,"m13*(mv03+3*mv02+3*mv01)"
   ,"m14*(mv04+4*mv03+6*mv02+4*mv01)"
   ,"m21*(mv11+mv01+mv10)"
   ,"m22*(mv12+mv02+2*mv11+2*mv01+mv10)"
   ,"m31*(mv21+2*mv11+mv01+mv20+2*mv10)"
   ,"m32*(mv22+2*mv12+mv02+2*mv21+4*mv11+2*mv01+mv20+2*mv10)"
   ,"m23*(mv13+mv03+3*mv12+3*mv02+3*mv11+3*mv01+mv10)"
   ,"m31*(mv31+3*mv21+3*mv11+mv01+mv30+3*mv20+3*mv10)")

   Mstar3 =
  c("m11*(mv10)"
   ,"m21*(mv20+2*mv10)"
   ,"m31*(mv30+3*mv20+3*mv10)"
   ,"m41*(mv40+4*mv30+6*mv20+4*mv10)"
   ,"m02*(mv01)"
   ,"m03*(mv02+2*mv01)"
   ,"m04*(mv03+3*mv02+3*mv01)"
   ,"m05*(mv04+4*mv03+6*mv02+4*mv01)"
   ,"m12*(mv11+mv01+mv10)"
   ,"m13*(mv12+mv02+2*mv11+2*mv01+mv10)"
   ,"m22*(mv21+2*mv11+mv01+mv20+2*mv10)"
   ,"m23*(mv22+2*mv12+mv02+2*mv21+4*mv11+2*mv01+mv20+2*mv10)"
   ,"m14*(mv13+mv03+3*mv12+3*mv02+3*mv11+3*mv01+mv10)"
   ,"m32*(mv31+3*mv21+3*mv11+mv01+mv30+3*mv20+3*mv10)")
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
 ##print(namess[whichnames])
 dims= rep('',14*2)
 for(i in 1:14)
 {
   for(j in whichnames)
   {
     if(MAT[i,j]!='')
     {
      dims[i] =paste0(dims[i],'+(',body(namess[j])[2],')*(',MAT[i,j],')')
      dims[i+14] =paste0(dims[i+14],'+(',body(namess[j])[2],')*(',MAT2[i,j],')')

     }
   }
 }
 ##print(MAT)
 #write.table(data.frame(dims),'OODEs.txt')

 #namess2 = c('Lam00','Lamy0','Nmu11','Nsig11','Nmu21','Nsig21','Nmu12','Nsig12','Nmu22','Nsig22','Lam10','Lamx2','Lamy1','Lamy2')
 namess2 = c('Lam00','Lam10','Lam01','Jmu1','Jmu2','Jsig11','Jsig12','Jsig22')
 checknames2 = rep(0,length(namess2))
 for(i in 1:length(namess2))
 {
    if(sum(objlist==namess2[i])==1){checknames2[i] = 1}
 }
 whichnames2 = which(checknames2==1)


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
     dims[i] = paste0(dims[i],'+(',body(namess2[1])[2],')*(',Mstar1[i],')')
  }
 }

 if(checknames2[2]==1)
 {
  for(i in 1:14)
  {
     dims[i] = paste0(dims[i],'+(',body(namess2[2])[2],')*(',Mstar2[i],')')
  }
 }

  if(checknames2[3]==1)
 {
  for(i in 1:14)
  {
     dims[i] = paste0(dims[i],'+(',body(namess2[3])[2],')*(',Mstar3[i],')')
  }
 }
 # if(checknames2[11]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[11])[2],')*(',Mstar1x[i],')')
 # }
 #}

 #  if(checknames2[12]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[12])[2],')*(',Mstar2x[i],')')
 # }
 #}
 #  if(checknames2[13]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[13])[2],')*(',Mstar1y[i],')')
 # }
 #}

 #  if(checknames2[14]==1)
 #{
 # for(i in 1:14)
 # {
 #    dims[i] = paste0(dims[i],'+(',body(namess2[14])[2],')*(',Mstar2y[i],')')
 # }
 #}

 #prem =
 #'
 #    mm1 = mu11
 #    mm2 = (mu11^2)+ (sig11^2)
 #    mm3 = (mu11^3)+ 3* (mu11^1)*(sig11^2)
 #    mm4 = (mu11^4)+ 6* (mu11^2)*(sig11^2)+3*(sig11^4)
 #
 #    oo1 = mu21
 #    oo2 = (mu21^2)+ (sig21^2)
 #    oo3 = (mu21^3)+ 3* (mu21^1)*(sig21^2)
 #    oo4 = (mu21^4)+ 6* (mu21^2)*(sig21^2)+3*(sig21^4)
 #
 #    nn1 = mu12
 #    nn2 = (mu12^2)+ (sig12^2)
 #    nn3 = (mu12^3)+ 3* (mu12^1)*(sig12^2)
 #    nn4 = (mu12^4)+ 6* (mu12^2)*(sig12^2)+3*(sig12^4)
 #
 #    pp1 = mu22
 #    pp2 = (mu22^2)+ (sig22^2)
 #    pp3 = (mu22^3)+ 3* (mu22^1)*(sig22^2)
 #    pp4 = (mu22^4)+ 6* (mu22^2)*(sig22^2)+3*(sig22^4)
 #'
  if(Jdist=='MVNormal')
 {
  prem =
 '
     mv10 = mu1                                                          ;
     mv20 = pow(mu1,2)+ sig11                                            ;
     mv30 = pow(mu1,3)+ 3* pow(mu1,1)*sig11                              ;
     mv40 = pow(mu1,4)+ 6* pow(mu1,2)*sig11+3*pow(sig11,2)               ;
     mv01 = mu2                                                          ;
     mv02 = pow(mu2,2)+ sig22                                            ;
     mv03 = pow(mu2,3)+ 3* pow(mu2,1)*sig22                              ;
     mv04 = pow(mu2,4)+ 6* pow(mu2,2)*sig22+3*pow(sig22,2)               ;
     mv11 = mu1*mu2+sig12                                                ;
     mv12 = mu1*pow(mu2,2)+2*mu2*sig12+mu1*sig22                         ;
     mv21 = pow(mu1,2)*mu2+2*mu1*sig12+mu2*sig11                         ;
     mv22 = pow(mu1,2)*pow(mu2,2)+pow(mu2,2)*sig11+pow(mu1,2)*sig22+4*mu1*mu2*sig12 +sig11*sig22+2*sig12*sig12;
     mv13 = mu1*pow(mu2,3)+3*pow(mu2,2)*sig12+3*mu1*mu2*sig22 +3*sig12*sig22;
     mv31 = mu2*pow(mu1,3)+3*pow(mu1,2)*sig12+3*mu1*mu2*sig11 +3*sig12*sig11;
 '
   }
     #premprem =rep('',4)
     #'Lam00','Lamy0','Nmu11','Nsig11','Nmu21','Nsig21','Nmu12','Nsig12','Nmu22','Nsig22','Lam10','Lamx2','Lamy1','Lamy2'
     #if(checknames2[3]==1){premprem[1]= paste0(' mu11=',body('Nmu11')[2],';',' sig11=',body('Nsig11')[2],';')}
     #if(checknames2[5]==1){premprem[3]= paste0(' mu21=',body('Nmu21')[2],';',' sig21=',body('Nsig21')[2],';')}
     #if(checknames2[7]==1){premprem[2]= paste0(' mu12=',body('Nmu12')[2],';',' sig12=',body('Nsig12')[2],';')}
     #if(checknames2[9]==1){premprem[4]= paste0(' mu22=',body('Nmu22')[2],';',' sig22=',body('Nsig22')[2],';')}
     #if(checknames2[3]==0){premprem[1]= paste0(' mu11=',0,';',' sig11=',0,';')}
     #if(checknames2[5]==0){premprem[3]= paste0(' mu21=',0,';',' sig21=',0,';')}
     #if(checknames2[7]==0){premprem[2]= paste0(' mu12=',0,';',' sig12=',0,';')}
     #if(checknames2[9]==0){premprem[4]= paste0(' mu22=',0,';',' sig22=',0,';')}

     premprem =rep('',4)
     if(Jdist=='MVNormal')
     {
        premprem[1]= paste0(' mu1=',body('Jmu1')[2],';\n',' mu2=',body('Jmu2')[2],';\n',' sig11=',body('Jsig11')[2],';\n',' sig12=',body('Jsig12')[2],';\n',' sig22=',body('Jsig22')[2],';\n')
        premprem[2:4]=''
     }


   prm=''
   for(i in 1:4)
   {
     prm = paste0(prm,premprem[i],'\n ')
   }

   if(checknames2[1]==1){b1 = paste0('(',body('Lam00')[2],')')}else{b1 = paste0('(',0,')')}
   if(checknames2[2]==1){b2 = paste0('(',body('Lam10')[2],')')}else{b2 = paste0('(',0,')')}
   if(checknames2[3]==1){b3 = paste0('(',body('Lam01')[2],')')}else{b3 = paste0('(',0,')')}
   #if(checknames2[12]==1){b4 = paste0('(',body('Lamx2')[2],')')}else{b4 = paste0('(',0,')')}
   #if(checknames2[13]==1){b5 = paste0('(',body('Lamy1')[2],')')}else{b5 = paste0('(',0,')')}
   #if(checknames2[14]==1){b6 = paste0('(',body('Lamy2')[2],')')}else{b6 = paste0('(',0,')')}
   b4 = paste0('(',0,')')
   b5 = paste0('(',0,')')
   b6 = paste0('(',0,')')

   odekernel=paste0('c(',dims[1],'\n,',dims[2],'\n,',dims[3],'\n,',dims[4],
                   '\n,',dims[5],'\n,',dims[6],'\n,',dims[7],'\n,',dims[8],'\n,',dims[9],'\n,',dims[10],
                   '\n,',dims[11],'\n,',dims[12],'\n,',dims[13],'\n,',dims[14],'\n,',dims[1+14],'\n,',dims[2+14],'\n,',dims[3+14],'\n,',dims[4+14],
                   '\n,',dims[5+14],'\n,',dims[6+14],'\n,',dims[7+14],'\n,',dims[8+14],'\n,',dims[9+14],'\n,',dims[10+14],
                   '\n,',dims[11+14],'\n,',dims[12+14],'\n,',dims[13+14],'\n,',dims[14+14],'\n,',b1,'\n,',b2,'\n,',b3,'\n,',b4,'\n,',b5,'\n,',b6,')')

    pr=
    '
    m00=1
    m10=x[1]
    m20=x[2]
    m30=x[3]
    m40=x[4]
    m01=x[5]
    m02=x[6]
    m03=x[7]
    m04=x[8]
    m11=x[9]
    m12=x[10]
    m21=x[11]
    m22=x[12]
    m13=x[13]
    m31=x[14]
    mm00=1
    mm10=x[1+14]
    mm20=x[2+14]
    mm30=x[3+14]
    mm40=x[4+14]
    mm01=x[5+14]
    mm02=x[6+14]
    mm03=x[7+14]
    mm04=x[8+14]
    mm11=x[9+14]
    mm12=x[10+14]
    mm21=x[11+14]
    mm22=x[12+14]
    mm13=x[13+14]
    mm31=x[14+14]


   k10 = m10
   k20 = m20-m10^2
   k30 = m30-3*m20*m10+2*m10^3
   k40 = m40 -4*m30*m10-3*m20^2+12*m20*m10^2-6*m10^4
   k01 = m01
   k02 = m02-  m01^2
   k03 = m03-3*m02*m01+2*m01^3
   k04 = m04-4*m03*m01-3*m02^2+12*m02*m01^2-6*m01^4
   k11 = m11-m10*m01
   k21 = m21-2*m11*m10-m20*m01+2*m10^2*m01
   k12 = m12-2*m11*m01-m02*m10+2*m01^2*m10
   k22 = m22-2*m21*m01-2*m12*m10-m20*m02-2*m11^2+8*m11*m01*m10+2*m02*m10^2+2*m20*m01^2-6*m10^2*m01^2
   k31 = m31-3*m21*m10-m30*m01-3*m20*m11+6*m11*m10^2+6*m20*m10*m01-6*m10^3*m01
   k13 = m13-3*m12*m01-m03*m10-3*m02*m11+6*m11*m01^2+6*m02*m01*m10-6*m01^3*m10

   m50 = 5*m10*m40+10*m20*m30-20*m10^2*m30-30*m10*m20^2+60*m10^3*m20-24*m10^5  #Correct
   m05 = 5*m01*m04+10*m02*m03-20*m01^2*m03-30*m01*m02^2+60*m01^3*m02-24*m01^5  #Correct
   m41 = +4*k31*k10+4*k30*k11+6*k20*k21+12*k20*k10*k11+6*k21*k10^2+4*k10^3*k11+m01*m40 #Correct
   m14 = +4*k13*k01+4*k03*k11+6*k02*k12+12*k02*k01*k11+6*k12*k01^2+4*k01^3*k11+m10*m04 #Correct
   m32 = +3*k20*k12+3*k21*k11+3*k21*k11+3*k22*k10+3*k10^2*k12+6*k11^2*k10+k01*k31+k02*k30+3*k02*k20*k10+3*k01*(k21*k10+k20*k11)+k02*k10^3+3*k01*k10^2+m31*m01
   m23 = +3*k02*k22+3*k12*k11+3*k12*k11+3*k22*k01+3*k01^2*k21+6*k11^2*k01+k10*k13+k20*k03+3*k20*k02*k01+3*k10*(k12*k01+k02*k11)+k20*k01^3+3*k10*k01^2+m13*m10

   kk10 = mm10
   kk20 = mm20-mm10^2
   kk30 = mm30-3*mm20*mm10+2*mm10^3
   kk40 = mm40 -4*mm30*mm10-3*mm20^2+12*mm20*mm10^2-6*mm10^4
   kk01 = mm01
   kk02 = mm02-  mm01^2
   kk03 = mm03-3*mm02*mm01+2*mm01^3
   kk04 = mm04-4*mm03*mm01-3*mm02^2+12*mm02*mm01^2-6*mm01^4
   kk11 = mm11-mm10*mm01
   kk21 = mm21-2*mm11*mm10-mm20*mm01+2*mm10^2*mm01
   kk12 = mm12-2*mm11*mm01-mm02*mm10+2*mm01^2*mm10
   kk22 = mm22-2*mm21*mm01-2*mm12*mm10-mm20*mm02-2*mm11^2+8*mm11*mm01*mm10+2*mm02*mm10^2+2*mm20*mm01^2-6*mm10^2*mm01^2
   kk31 = mm31-3*mm21*mm10-mm30*mm01-3*mm20*mm11+6*mm11*mm10^2+6*mm20*mm10*mm01-6*mm10^3*mm01
   kk13 = mm13-3*mm12*mm01-mm03*mm10-3*mm02*mm11+6*mm11*mm01^2+6*mm02*mm01*mm10-6*mm01^3*mm10

   mm50 = 5*mm10*mm40+10*mm20*mm30-20*mm10^2*mm30-30*mm10*mm20^2+60*mm10^3*mm20-24*mm10^5  #Correct
   mm05 = 5*mm01*mm04+10*mm02*mm03-20*mm01^2*mm03-30*mm01*mm02^2+60*mm01^3*mm02-24*mm01^5  #Correct
   mm41 = +4*kk31*kk10+4*kk30*kk11+6*kk20*kk21+12*kk20*kk10*kk11+6*kk21*kk10^2+4*kk10^3*kk11+mm01*mm40 #Correct
   mm14 = +4*kk13*kk01+4*kk03*kk11+6*kk02*kk12+12*kk02*kk01*kk11+6*kk12*kk01^2+4*kk01^3*kk11+mm10*mm04 #Correct
   mm32 = +3*kk20*kk12+3*kk21*kk11+3*kk21*kk11+3*kk22*kk10+3*kk10^2*kk12+6*kk11^2*kk10+kk01*kk31+kk02*kk30+3*kk02*kk20*kk10+3*kk01*(kk21*kk10+kk20*kk11)+kk02*kk10^3+3*kk01*kk10^2+mm31*mm01
   mm23 = +3*kk02*kk22+3*kk12*kk11+3*kk12*kk11+3*kk22*kk01+3*kk01^2*kk21+6*kk11^2*kk01+kk10*kk13+kk20*kk03+3*kk20*kk02*kk01+3*kk10*(kk12*kk01+kk02*kk11)+kk20*kk01^3+3*kk10*kk01^2+mm13*mm10




    '
   odekernel=paste0('{ \n',pr,prm,prem,odekernel,'}')
  # write(odekernel,'res.txt')
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

   type.sol ="                  GENERALIZED QUADRATIC DIFFUSON"
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
       if(sum(objlist==namess4[i]))
       {
         indnames[i]=TRUE
         namess4[i]=paste0(namess4[i],' : ',trim(body(namess4[i])[2]))
        }
   }
   indnames2 = rep(0,5)
   for(j in whichnames2)
   {
       if(sum(objlist==namess2[j]))
       {
         indnames2[j]=TRUE
         namess2[j]=paste0(namess2[j],' : ',trim(body(namess2[j])[2]))
       }
   }
   #prior.list = trim(prior.list)
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
         namess2[c(4:8)][which(indnames2[4:8]==T)])
   Info=data.frame(matrix(Info,length(Info),1))
   colnames(Info)=''
   if(print.output)
   {
     print(Info,row.names = FALSE,right=F)
   }

       N=round((t-s)/delt)
       ttt=seq(s,t,by=delt)
       seq1=c(1,2,3,4,0,0,0,0,1,1,2,2,1,3)
       seq2=c(0,0,0,0,1,2,3,4,1,2,1,2,3,1)
       x00= rep(0,14*2+6)
       for(i in 1:length(seq1))
       {
          x00[i] = (Xs^seq1[i])*(Ys^seq2[i])
          x00[i+14] = (Xs^seq1[i])*(Ys^seq2[i])
       }
       pb <- txtProgressBar(1,2*N,1,style = 1,width = 65)
       #print(x00)
       solver=function(x00)
       {
         MM=matrix(0,2*14+6,N+1)
         MA= x00
         tt=s

         MM[,1]=MA
         ii=1
         while(ii < N+1)
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

     MM=solver(x00)
     #print(length(ff(x00,tt)))
     #print(ff(x00,tt))
     #return(0)
    # windows()
    # plot(MM[1,],type='n',ylim=range(log(MM)))
    # for(i in 1:14)
    # {
    #     lines(log(MM[i,]),col=2)
    #     lines(log(MM[i+14,]),col=4)
    # }
    nnn = length(Xt)



     EE =MM[1:14,]*0
     for(i in 1:14)
     {
        EE[i,] = (MM[i,]-MM[i+14,]*exp(-MM[29,]-MM[30,]-MM[31,]-MM[32,]-MM[33,]-MM[34,]))/(1-exp(-MM[29,]-MM[30,]-MM[31,]-MM[32,]-MM[33,]-MM[34,]))
     }

  if(eval.density)
  {
     if(Dtype=='Saddlepoint')
  {
    n1=dim(MM)[1]
    n2=dim(MM)[2]

    kkk=0
    DDD=array(0,dim=c(nnn,nnn,n2))

    a=rep(0,nnn*nnn)
    b=rep(0,nnn*nnn)

    saddlepoint=function(X1,X2,m)
    {
          m00=1
          m10=m[1]
          m20=m[2]
          m30=m[3]
          m40=m[4]
          m01=m[5]
          m02=m[6]
          m03=m[7]
          m04=m[8]
          m11=m[9]
         m12=m[10]
         m21=m[11]
         m22=m[12]
         m13=m[13]
         m31=m[14]

         k10 = m10
         k20 = m20-m10^2
         k30 = m30-3*m20*m10+2*m10^3
         k40 = m40 -4*m30*m10-3*m20^2+12*m20*m10^2-6*m10^4
         k01 = m01
         k02 = m02-  m01^2
         k03 = m03-3*m02*m01+2*m01^3
         k04 = m04-4*m03*m01-3*m02^2+12*m02*m01^2-6*m01^4
         k11 = m11-m10*m01
         k21 = m21-2*m11*m10-m20*m01+2*m10^2*m01
         k12 = m12-2*m11*m01-m02*m10+2*m01^2*m10
         k22 = m22-2*m21*m01-2*m12*m10-m20*m02-2*m11^2+8*m11*m01*m10+2*m02*m10^2+2*m20*m01^2-6*m10^2*m01^2
         k31 = m31-3*m21*m10-m30*m01-3*m20*m11+6*m11*m10^2+6*m20*m10*m01-6*m10^3*m01
         k13 = m13-3*m12*m01-m03*m10-3*m02*m11+6*m11*m01^2+6*m02*m01*m10-6*m01^3*m10

         ffab=function(xm1,xm2,a=rep(0,nnn*nnn),b=rep(0,nnn*nnn))
        {
           XMAT1=xm1
          XMAT2=xm2


        k=0
        abser=rep(0.1,nnn*nnn)
        while((k<5000)&(sum(abser>0.005)>0))
        {
          # K= k10*a+k01*b+1/2*k20*a^2+1/2*k02*b^2+1/6*k30*a^3+1/6*k03*b^3+1/24*k40*a^4+1/24*k04*b^4+k11*a*b+1/2*k12*a*b^2+1/2*k21*a^2*b+1/6*a*b^3*k13+1/6*a^3*b*k31+1/4*a^2*b^2*k22
          gg=k10+k20*a+1/2*k30*a^2+1/6*k40*a^3 +k11*b +1/2*k12*b^2+k21*a*b+1/6*b^3*k13+1/2*a^2*b*k31+1/2*a*b^2*k22-XMAT1
          hh=k01+k02*b+1/2*k03*b^2+1/6*k04*b^3 +k11*a +k12*a*b+1/2*k21*a^2+1/2*a*b^2*k13+1/6*a^3*k31+1/2*a^2*b*k22-XMAT2
          gg1=k20+k30*a+1/2*k40*a^2+k21*b+a*b*k31+1/2*b^2*k22
          gg2=k11 +k12*b+k21*a+1/2*b^2*k13+1/2*a^2*k31+a*b*k22
          hh1=k11 +k12*b+k21*a+1/2*b^2*k13+1/2*a^2*k31+a*b*k22
          hh2=k02+k03*b+1/2*k04*b^2 +k12*a+a*b*k13+1/2*a^2*k22
          DK2=gg1*hh2-gg2*hh1
          ar=(hh*gg2-gg*hh2)/(gg1*hh2-gg2*hh1)
          br=(hh*gg1-gg*hh1)/(gg1*hh2-gg2*hh1)
          anew  =a+ar
          bnew  =b-br
          abser =((anew-a)^2+(bnew-b)^2)
          itrue =(abser>0.001)
          ifalse=(abser<=0.001)
          #print(itrue)
          a=anew#*ifalse+anew*itrue
          b=bnew#*ifalse+bnew*itrue
          k=k+1
        }

        return(list(a=a,b=b,k=k))
      }

      fff=function(a,b,xmat1,xmat2)
      {
        K= k10*a+k01*b+1/2*k20*a^2+1/2*k02*b^2+1/6*k30*a^3+1/6*k03*b^3+1/24*k40*a^4+1/24*k04*b^4+k11*a*b+1/2*k12*a*b^2+1/2*k21*a^2*b+1/6*a*b^3*k13+1/6*a^3*b*k31+1/4*a^2*b^2*k22
        gg=k10+k20*a+1/2*k30*a^2+1/6*k40*a^3 +k11*b +1/2*k12*b^2+k21*a*b+1/6*b^3*k13+1/2*a^2*b*k31+1/2*a*b^2*k22-xmat1
        hh=k01+k02*b+1/2*k03*b^2+1/6*k04*b^3 +k11*a +k12*a*b+1/2*k21*a^2+1/2*a*b^2*k13+1/6*a^3*k31+1/2*a^2*b*k22-xmat2
        gg1=k20+k30*a+1/2*k40*a^2+k21*b+a*b*k31+1/2*b^2*k22
        gg2=k11 +k12*b+k21*a+1/2*b^2*k13+1/2*a^2*k31+a*b*k22
        hh1=k11 +k12*b+k21*a+1/2*b^2*k13+1/2*a^2*k31+a*b*k22
        hh2=k02+k03*b+1/2*k04*b^2 +k12*a+a*b*k13+1/2*a^2*k22
        DK2=gg1*hh2-gg2*hh1
        return(exp(K-a*xmat1-b*xmat2)/(2*pi)/sqrt(abs(DK2)))
      }
             dett = (k02*k20-k11^2)
      a    = (k20*(k10-X1)-k11*(k01-X2))/dett
      b    = (-k11*(k10-X1)+k02*(k01-X2))/dett
      cc=ffab(X1,X2,a,b)
      return(t(matrix(fff(cc$a,cc$b,X1,X2),nnn,nnn,byrow=T)))
     }




    for(i in 10:n2)
    {
      setTxtProgressBar(pb,i+n2," "," ")
      kkk=kkk+1
      x=MM[,i]

      p=exp(-MM[29,i]-MM[30,i]-MM[31,i]-MM[32,i]-MM[33,i]-MM[34,i])
      #print(p)
      #print(EE[,i])
      DDD[,,i] =  saddlepoint(X1,X2,MM[1:14+14,i])*p+saddlepoint(X1,X2,EE[,i])*(1-p)
   }
   }

   if(Dtype=='Edgeworth')
   {

        ffff=function(xx,yy,m)
       {
        X=rep(xx,length(yy))
        Y=rep(yy,each=length(xx))

         H=function(x,i)
         {
            switch(i+1,
            {1},
            {x},
            {x^2-1},
            {x^3-3*x},
            {x^4-6*x^2+3},
            {x^5-10*x^3+15*x},
            {x^6-15*x^4+45*x^2-15})
         }

          m00=1
          m10=m[1]
          m20=m[2]
          m30=m[3]
          m40=m[4]
          m01=m[5]
          m02=m[6]
          m03=m[7]
          m04=m[8]
          m11=m[9]
         m12=m[10]
         m21=m[11]
         m22=m[12]
         m13=m[13]
         m31=m[14]

         k10 = m10
         k20 = m20-m10^2
         k30 = m30-3*m20*m10+2*m10^3
         k40 = m40 -4*m30*m10-3*m20^2+12*m20*m10^2-6*m10^4
         k01 = m01
         k02 = m02-  m01^2
         k03 = m03-3*m02*m01+2*m01^3
         k04 = m04-4*m03*m01-3*m02^2+12*m02*m01^2-6*m01^4
         k11 = m11-m10*m01
         k21 = m21-2*m11*m10-m20*m01+2*m10^2*m01
         k12 = m12-2*m11*m01-m02*m10+2*m01^2*m10
         k22 = m22-2*m21*m01-2*m12*m10-m20*m02-2*m11^2+8*m11*m01*m10+2*m02*m10^2+2*m20*m01^2-6*m10^2*m01^2
         k31 = m31-3*m21*m10-m30*m01-3*m20*m11+6*m11*m10^2+6*m20*m10*m01-6*m10^3*m01
         k13 = m13-3*m12*m01-m03*m10-3*m02*m11+6*m11*m01^2+6*m02*m01*m10-6*m01^3*m10

          x = (X-k10)/sqrt(k20)
          y = (Y-k01)/sqrt(k02)
          rho20 =1
          rho30 =k30/(sqrt(k20)^3)
          rho40 =k40/(sqrt(k20)^4)

          rho02 =1
          rho03 =k03/(sqrt(k02)^3)
          rho04 =k04/(sqrt(k02)^4)

          rho11 =k11/(sqrt(k20)*sqrt(k02))
          rho12 =k12/(sqrt(k20)*sqrt(k02)^2)
          rho21 =k21/(sqrt(k20)^2*sqrt(k02))
          rho22 =k22/(sqrt(k20)^2*sqrt(k02)^2)
          rho13 =k13/(sqrt(k20)^1*sqrt(k02)^3)
          rho31 =k31/(sqrt(k20)^3*sqrt(k02)^1)


          rho = k11/sqrt(k02*k20)
          dmvnorm <- function(x,y)
          {
          	 1/(2*pi*sqrt(k20*k02)*sqrt(1 - rho^2))*exp(-((x-k10)^2/k20 -2*rho*(x-k10)*(y-k01)/sqrt(k20 * k02)+(y-k01)^2/k02)/(2*(1-rho^2)))
          }

        A=H(x,3)*rho30+3*H(x,2)*H(y,1)*rho21+3*H(x,1)*H(y,2)*rho12+H(y,3)*rho03
        B=H(x,4)*rho40+4*H(x,3)*H(y,1)*rho31+6*H(x,2)*H(y,2)*rho22+4*H(x,1)*H(y,3)*rho13+H(y,4)*rho04
        CC=H(x,6)*rho30^2+6*H(x,5)*H(y,1)*rho21*rho30+6*H(x,4)*H(y,2)*rho12*rho30+2*H(x,3)*H(y,3)*rho03*rho30+9*H(x,4)*H(y,2)*rho21^2+18*H(x,3)*H(y,3)*rho12*rho21+6*H(x,2)*H(y,4)*rho03*rho21+9*H(x,2)*H(y,4)*rho12^2+6*H(x,1)*H(y,5)*rho03*rho12+H(y,6)*rho03^2
        val=pmax(dmvnorm(X,Y)*(1+1/6*A+1/24*B+1/72*CC),0)


       return(val)
       }

       DDD=array(0,dim=c(nnn,nnn,dim(MM)[2]))
       DDD[1,1,1] = 1
       for(i in 2:dim(MM)[2])
       {
         #print(exp(-MM[29,i]-MM[30,i]-MM[31,i]-MM[32,i]-MM[33,i]-MM[34,i]))
         DDD[,,i]=exp(-MM[29,i]-MM[30,i]-MM[31,i]-MM[32,i]-MM[33,i]-MM[34,i])*ffff(Xt,Yt,MM[15:28,i])+(1-exp(-MM[29,i]-MM[30,i]-MM[31,i]-MM[32,i]-MM[33,i]-MM[34,i]))*ffff(Xt,Yt,EE[,i])
         setTxtProgressBar(pb,i+N," "," ")
       }
   }
   close(pb)
    DD1=matrix(0,nnn,dim(MM)[2])
           saddlep=function(xx,m)
       {
        k =m*0
        k[1,]=                                          m[1,]
        k[2,]=                               m[2,]-1*k[1,]*m[1,]
        k[3,]=                m[3,]-1*k[1,]*m[2,]-2*k[2,]*m[1,]
        k[4,]= m[4,]-1*k[1,]*m[3,]-3*k[2,]*m[2,]-3*k[3,]*m[1,]

        p=1/3*(3*(k[4,]/6)*k[2,] - ((k[3,]/2)^2))/((k[4,]/6)^2)
        q=1/27*(27*((k[4,]/6)^2)*(k[1,]-xx) - 9*(k[4,]/6)*(k[3,]/2)*k[2,] + 2*((k[3,]/2)^3))/((k[4,]/6)^3)
        chk=(q^2)/4 + (p^3)/27
        th=-(k[3,]/2)/(3*(k[4,]/6))+(-q/2+sqrt(chk))^(1/3)-(q/2+sqrt(chk))^(1/3)

        K=k[1,]*th+(k[2,]*th^2)/2+(k[3,]*th^3)/6 +(k[4,]*th^4)/24
        K1=k[1,]+(k[2,]*th)+(k[3,]*th^2)/2+(k[4,]*th^3)/6
        K2=k[2,]+(k[3,]*th)+(k[4,]*th^2)/2
        K3=k[3,]+(k[4,]*th)
        K4=k[4,]
        return(1/sqrt(2*pi*(K2))*exp(K-th*K1))
       }
      for(i in 1:length(Xt))
    {

      DD1[i,]=saddlep(Xt[i],MM[1:4,])
    }
    DD2=matrix(0,nnn,dim(MM)[2])
      for(i in 1:length(Yt))
    {
     DD2[i,]=saddlep(Yt[i],MM[1:4+4,])
    }
   }
   if(!eval.density)
   {
     DDD=NULL;DD1=NULL;DD2=NULL;setTxtProgressBar(pb,2*N," "," ");
   }
   rownames(MM) =c("m10","m20","m30","m40","m01","m02","m03","m04","m11","m12","m21","m22","m13","m31","mm10","mm20","mm30","mm40","mm01","mm02","mm03","mm04","mm11","mm12","mm21","mm22","mm13","mm31",'ILam00','ILam10','ILam01','ILamx2','ILamy1','ILamy2')

   retlist=list(density=DDD,time=seq(s,t,delt),Xt=Xt,Yt=Yt,moments=MM[1:31,],Xmarginal=DD1,Ymarginal=DD2,excess_moments=EE,zero_jump_prob=exp(-MM[29,]-MM[30,]-MM[31,]-MM[32,]-MM[33,]-MM[34,]))
   class(retlist) = "BiJGQD.density"
   return(retlist)
 }







