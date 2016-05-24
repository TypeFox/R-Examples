`strikeslip.fault` <-
function(anim= seq(from=0, to=1, by=.1), KAPPA = 2,  Light=c(45,45) )
  {
    if(missing(anim)) {  anim= seq(from=0, to=1, by=.1)  }

    if(missing(KAPPA)) {    KAPPA = 2 }
    if(missing(Light)) { Light=c(45, 45) }

    
    block1 = matrix(c(0,0,0,
      1,0,0,
      1,0.5,0,
      0,0.5,0,
      0,0,-2,
      1,0,-2,
      1,0.5,-2,
      0,0.5,-2), byrow=TRUE, ncol=3)

    Bblock1 = makeblock3D(block1)

    i = max(anim)
    angx = -40
    angz = 20
    alpha = angz
    gamma = angx
    beta  = 0
    R3 = ROTX(gamma)
    R2 = ROTY(beta)
    R1 = ROTZ(alpha)
    T =  TRANmat(-1, 0, 0 )
    M =     R1  %*% R2  %*%  R3  %*% T  

    T2 =  TRANmat(-i, 0.5, 0 )
    MT =       T2 %*%   R1  %*% R2  %*%   R3 %*% T   

    Z1 =  PROJ3D(Bblock1$aglyph, M=MT,  anorms=Bblock1$anorm , zee=c(0,0,1))

    Z2 =     PROJ3D(Bblock1$aglyph, M=M,  anorms=Bblock1$anorm , zee=c(0,0,1))
#######################         #######################

    
    angx = -40
    angz = 20-90
    alpha = angz
    gamma = angx
    beta  = 0
    R3 = ROTX(gamma)
    R2 = ROTY(beta)
    R1 = ROTZ(alpha)
    T =  TRANmat(1, 0, 0 )
    M =      R1  %*% R2  %*%  R3 %*% T
    T2 =  TRANmat(-i, 0.5, 0 )
    MT =         T2 %*%   R1  %*% R2  %*%  R3   %*% T

    Z3 =  PROJ3D(Bblock1$aglyph, M=MT,  anorms=Bblock1$anorm , zee=c(0,0,1))

    Z4 =     PROJ3D(Bblock1$aglyph, M=M,  anorms=Bblock1$anorm , zee=c(0,0,1))
#######################         #######################

    RangesX = range(c(attr(Z1, "RangesX"), attr(Z2, "RangesX"), attr(Z3, "RangesX") , attr(Z4, "RangesX")))

    RangesY = range(c(attr(Z1, "RangesY"), attr(Z2, "RangesY"), attr(Z3, "RangesY") , attr(Z4, "RangesY")))
#######################

    TFanim = TRUE 
    while(TFanim)
      {
        

        for(i in  anim  )
          {
            dkm = i
            dkm2 = -i
            
            angx = -40
            angz = 20
            alpha = angz
            gamma = angx
            beta  = 0
            R3 = ROTX(gamma)
            R2 = ROTY(beta)
            R1 = ROTZ(alpha)
            T =  TRANmat(-1, 0, 0 )
            M =     R1  %*% R2  %*%  R3  %*% T  

####    plot(c(-2,2) , c(-2,2) , type='n', asp=1, ann=FALSE, axes=FALSE)
            plot( RangesX , RangesY,  type='n', asp=1, ann=FALSE, axes=FALSE)
            T2 =  TRANmat(dkm2, 0.5, 0 )
            MT =       T2 %*%   R1  %*% R2  %*%   R3 %*% T   

            phong3D(Bblock1$aglyph, M=MT,  anorms=Bblock1$anorm , zee=c(0,0,1), Light=Light, col=rgb(1,.7,.7), border="black")

            phong3D(Bblock1$aglyph, M=M,   anorms=Bblock1$anorm , zee=c(0,0,1), Light=Light, col="white", border="black")
            ## abline(v=0, h=0)

            angx = -40
            angz = 20-90
            alpha = angz
            gamma = angx
            beta  = 0
            R3 = ROTX(gamma)
            R2 = ROTY(beta)
            R1 = ROTZ(alpha)
            T =  TRANmat(1, 0, 0 )
            M =      R1  %*% R2  %*%  R3 %*% T
            T2 =  TRANmat(dkm2, 0.5, 0 )
            MT =         T2 %*%   R1  %*% R2  %*%  R3   %*% T

            phong3D(Bblock1$aglyph, M=MT,  anorms=Bblock1$anorm , zee=c(0,0,1), Light=Light, col=rgb(1,.7,.7), border="black")

            phong3D(Bblock1$aglyph, M=M,  anorms=Bblock1$anorm , zee=c(0,0,1), Light=Light, col="white", border="black")


            Sys.sleep(0.1)
            ## locator(1)
          }
        if(length(anim)>1) { TFanim = TRUE  } else { TFanim = FALSE }
      }



  }

