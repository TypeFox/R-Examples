`normal.fault` <-
function(ANG=(45), anim=seq(from=0, to=1, by=.1) , KAPPA = 4,  Light=c(45,45) )
{
###########   usage:   normal.fault()

######

  if(missing(ANG)) { ANG= (45) }
  if(missing(anim)) {  anim= seq(from=0, to=1, by=.1)  }

     if(missing(KAPPA)) {  KAPPA = 4 }
  if(missing(Light)) { Light=c(-20, 20) }

  DEGRAD = pi/180
  y1 = 1.5
  
  y2 = y1 - 1/tan((ANG)*DEGRAD)

  
  z1 = 1
  x1 = 1

  
  Ablock1 = matrix(c(0,0,0,
    1,0,0,
    1,y1,0,
    0,y1,0,
    0,0,-1,
    1,0,-1,
    1,y2,-1,
    0,y2,-1), byrow=TRUE, ncol=3)
  rblock1 = cbind(Ablock1, rep(1, length(Ablock1[,1])))

  AR2 = ROTX(-180)

  ##   print(c(y1, y2) )

  

  Ty =  TRANmat(0,y1+y2,-1 )

  ##  Ty =  TRANmat(0,0,0 )


  Ablock2 = rblock1 %*% AR2  %*% Ty  

  Ablock2  = Ablock2[,1:3]
  Ablock2[abs(Ablock2)<1.224647e-10] = 0.0

  ## block2[,2] = block2[,2] 
  
  if(ANG>90)
    {
      Nblock2 = makeblock3D(Ablock1)
      Nblock1 = makeblock3D(Ablock2)
    }else{
      Nblock1 = makeblock3D(Ablock1)
      Nblock2 = makeblock3D(Ablock2)
  

    }


  if(ANG<90){
  ANG2 = ANG+90
} else {
     ANG2 = ANG-90

   }


  
    y1 = 1.5
  
  y2 = y1 - 1/tan((ANG2)*DEGRAD)

  
  z1 = 1
  x1 = 1

  
  Bblock1 = matrix(c(0,0,0,
    1,0,0,
    1,y1,0,
    0,y1,0,
    0,0,-1,
    1,0,-1,
    1,y2,-1,
    0,y2,-1), byrow=TRUE, ncol=3)
  
  rblock1 = cbind(Bblock1, rep(1, length(Bblock1[,1])))

  AR2 = ROTX(-180)

  Ty =  TRANmat(0,y1+y2,-1 )

  ##  Ty =  TRANmat(0,0,0 )


  Bblock2 = rblock1 %*% AR2  %*% Ty  

  Bblock2  = Bblock2[,1:3]
  Bblock2[abs(Bblock2)<1.224647e-10] = 0.0

  ## block2[,2] = block2[,2] 
  
  if(ANG2>90)
    {
      Oblock2 = makeblock3D(Bblock1)
      Oblock1 = makeblock3D(Bblock2)
    } else
    {
      Oblock1 = makeblock3D(Bblock1)
      Oblock2 = makeblock3D(Bblock2)
  

    }

  angx = -45
  angz = -45

  beta  = 0
  R3 = ROTX(angx)
  R1 = ROTZ(angz)
  
  Tshift1 =  TRANmat(-3,0, 0 )
  Tshift2 =  TRANmat(1,0, 0 )


########  SANG = atan2( (Ablock2[4,3]-Ablock2[5,3])  ,   (Ablock2[4,3]-Ablock2[5,3])  )
  SANG = (ANG)*DEGRAD

  shiftx =  sin(0)*sin(SANG)
  shifty =  cos(0)*sin(SANG)
  shiftz =  cos(SANG)

  SANG2 = (ANG2)*DEGRAD

  shiftx2 =  sin(0)*sin(SANG2)
  shifty2 =  cos(0)*sin(SANG2)
  shiftz2 =  cos(SANG2)


###########   set up scene:
  dkm = max(anim)
          dkm2 = -dkm
          T2 =   TRANmat(dkm*shiftx, dkm*shifty, dkm*shiftz ) 
          T1 =  TRANmat(0,0,0 ) 
          MN =    T2 %*%  R1 %*%  R3    %*% Tshift1
          M =    T1 %*%  R1  %*%  R3   %*% Tshift1
          TO2 =   TRANmat(dkm2*shiftx2, dkm2*shifty2, dkm2*shiftz2 ) 
          TO1 =  TRANmat(0,0,0 ) 
          MNO =    TO2 %*%  R1 %*%  R3    %*% Tshift2
          MO =    TO1 %*%  R1  %*%  R3   %*% Tshift2

  Z1 = PROJ3D(Nblock2$aglyph, M=MN,  anorms=Nblock2$anorm , zee=c(0,0,1))
  Z2 = PROJ3D(Nblock1$aglyph, M=M,  anorms=Nblock1$anorm ,  zee=c(0,0,1))
  Z3 = PROJ3D(Oblock2$aglyph, M=MNO,  anorms=Oblock2$anorm ,  zee=c(0,0,1))
  Z4 = PROJ3D(Oblock1$aglyph, M=MO,  anorms=Oblock1$anorm ,  zee=c(0,0,1))


  RangesX = range(c(attr(Z1, "RangesX"), attr(Z2, "RangesX"), attr(Z3, "RangesX") , attr(Z4, "RangesX")))

  RangesY = range(c(attr(Z1, "RangesY"), attr(Z2, "RangesY"), attr(Z3, "RangesY") , attr(Z4, "RangesY")))

       
  TFanim = TRUE 
  while(TFanim)
    {
      for(i in anim )
        {
          dkm = i
          dkm2 = -i
           
          
###
###T2 =   TRANmat(0, -dkm, 0 ) 
          T2 =   TRANmat(dkm*shiftx, dkm*shifty, dkm*shiftz ) 

          T1 =  TRANmat(0,0,0 ) 
###  print(dkm)

          MN =    T2 %*%  R1 %*%  R3    %*% Tshift1
          M =    T1 %*%  R1  %*%  R3   %*% Tshift1



          TO2 =   TRANmat(dkm2*shiftx2, dkm2*shifty2, dkm2*shiftz2 ) 

          TO1 =  TRANmat(0,0,0 ) 
###  print(dkm)

          MNO =    TO2 %*%  R1 %*%  R3    %*% Tshift2
          MO =    TO1 %*%  R1  %*%  R3   %*% Tshift2

           plot( RangesX, RangesY, type='n', asp=1, ann=FALSE, axes=FALSE)
     ###      plot(c(-4,4) , c(-4,4) , type='n', asp=1, ann=FALSE, axes=FALSE)
###  abline(v=0,h=0, col=grey(0.85) , lty=2)

          phong3D(Nblock2$aglyph, M=MN,  anorms=Nblock2$anorm , Light = Light, zee=c(0,0,1), col=rgb(1,.7,.7) , border="black")

          phong3D(Nblock1$aglyph, M=M,  anorms=Nblock1$anorm , Light = Light, zee=c(0,0,1), col=grey(0.9) , border="black")


           phong3D(Oblock2$aglyph, M=MNO,  anorms=Oblock2$anorm , Light = Light, zee=c(0,0,1), col=rgb(1,.7,.7) , border="black")

          phong3D(Oblock1$aglyph, M=MO,  anorms=Oblock1$anorm , Light = Light, zee=c(0,0,1), col=grey(0.9) , border="black")


           

###X1 = cbind(0,1,0,1)  %*% M
###X2 = cbind(0,c(1:5),0,1)  %*% MN
###points(X1[,1], X1[,2] , pch=2, col='green')
###points(X2[,1], X2[,2] , pch=3, col='blue')

### locator(1)

          Sys.sleep(.1)
          
        }
      if(length(anim)>1) { TFanim = TRUE  } else { TFanim = FALSE }
      
    }


}

