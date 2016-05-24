`zoeppritz` <-
function(icoef, vp1, vp2, vs1, vs2, rho1, rho2, incw)
  {
    ##X##  calculate the zoeppritz equations
    ##X##    #   % calculate zoeppritz according to young and braille
    ##X## #  zoeppritz(icoef, 3.0, 6.0, , vs2, rho1, rho2, incw)
    ##X## #  zoeppritz(icoef, alpha1, alpha2, beta1, beta2, rho1,rho2,   incw)
    ##X## # alpha1 = 3.0;
    ##X## # beta1 = 1.7321;
    ##X## # rho1 = 2.69;

    ##X## # alpha2 = 6.0;
    ##X## # beta2 = 3.46;
    ##X## # rho2 = 2.91;
    ##X## #   incw=2;    icoef=1
    ##X##  #

    ##X## #    alpha1 = 3.0; beta1 = 1.7321;rho1 = 2.69; alpha2 = 6.0; beta2 = 3.46; rho2 = 2.91;

    ##X## #  icoef=1; incw=2; vp1=alpha1;  vp2=alpha2;  vs1=beta1; vs2=beta2;  rho1=rho1;  rho2=rho2
    
    ##X## #  zoeppritz(icoef, alpha1, alpha2, beta1, beta2, rho1,rho2,   incw)
    
    ##X## #   plotzpr( "Amplitude" , alpha1, alpha2, beta1, beta2, rho1 ,rho2, "SV", "S");


    if(incw==1)
      {
        vel = vp1;
      }
    if(incw==2)
      {
        vel = vs1;
      }

    angle = seq(from=0.1, to=89.7, length=90)
    alpha = angle*pi/180;

    fp = matrix(rep(0,4),  ncol=1,nrow=4);
    K = 4* length(alpha)
    
    r = matrix( rep(0,K), ncol=4, nrow=length(alpha));
    rmat = matrix( rep(0,K), ncol=4, nrow=length(alpha));
    
    fe = matrix( rep(0,K), ncol=4, nrow=length(alpha));

    rra = matrix( rep(0,K), ncol=4, nrow=length(alpha));
    ria = matrix( rep(0,K), ncol=4, nrow=length(alpha));
    ang = matrix( rep(0,K), ncol=4, nrow=length(alpha));


    theta  = sin(alpha)/vel;

    fac1 = rep(1,length(alpha));
    fac2 =rep(1,length(alpha));
    fac3 = rep(1,length(alpha));
    fac4 = rep(1,length(alpha));

    thetasq = theta*theta;
    qa = 2.0 *(rho2*vs2*vs2 - rho1*vs1*vs1);
    t1 = vp1*vp1*thetasq;
    t2 = vs1*vs1*thetasq;
    t3 = vp2*vp2*thetasq;
    t4 = vs2*vs2*thetasq;

                                        # %%%%%%%%%%%%%%%%%%%%    test for critical P refl 
    b1 = rep(0, length(alpha));
    a1 =rep(0, length(alpha));

    flag =  	 theta>1.0/vp1  ;

    b1[flag] = -sqrt(t1[flag]-1.0);
    a1[flag] = 0.0;
    fac1[flag] = 0.0;

    a1[!flag]=sqrt(1.0-t1[!flag]);
    b1[!flag] = 0.0;

                                        #   %%%%%%%%%%%%%%%%%%%%    test for critical S refl 
    b2 =  rep(0, length(alpha));
    a2 =  rep(0, length(alpha));

    flag =  	 theta>1.0/vs1  ;

    b2[flag] = -sqrt(t2[flag]-1.0);
    a2[flag] = 0.0;
    fac2[flag] = 0.0;

    a2[!flag]=sqrt(1.0-t2[!flag]);
    b2[!flag] = 0.0;
                                        #   %%%%%%%%%%%%%%%%%%%%    test for critical  p refraction 

    b3 =  rep(0, length(alpha));
    a3 = rep(0, length(alpha));
    flag = theta>1.0/vp2  ;

    b3[flag] = -sqrt(t3[flag]-1.0);
    a3[flag] = 0.0;
    fac3[flag] = 0.0;

    a3[!flag]=sqrt(1.0-t3[!flag]);
    b3[!flag] = 0.0;


                                        #   %%%%%%%%%%%%%%%%%%%%    test for critical  s refraction 
    b4 =  rep(0, length(alpha));
    a4 =  rep(0, length(alpha));
    flag =  theta>1.0/vs2  ;

    b4[flag] = -sqrt(t4[flag]-1.0);
    a4[flag] = 0.0;
    fac4[flag] = 0.0;

    a4[!flag]=sqrt(1.0-t4[!flag]);
    b4[!flag] = 0.0;

                                        #   %%%%%%%%%%%%%%%%%%%%

    x = rho2-(qa*thetasq);
    y = rho1+(qa*thetasq);
    z = rho2-rho1-(qa*thetasq) ;

    p1 = complex(real=a1, imaginary=b1);
    p2 = complex(real=a2, imaginary=b2);
    p3 = complex(real=a3, imaginary=b3);
    p4 = complex(real=a4, imaginary=b4);



    d = vp1*vp2*vs1*vs2*thetasq*z*z+ vp2*vs2*p1*p2*x*x+ vp1*vs1*p3*p4*y*y+
      rho1*rho2*(vs1*vp2*p1*p4+vp1*vs2*p2*p3)+  qa*qa*thetasq*p1*p2*p3*p4 ;

    fp = matrix(rep(1,4), ncol=2, nrow=2);

    if(incw==1)
      {
        
                                        #case 1

        r[,1 ] = -1.0 +2.0*p1*(vp2*vs2*p2*x*x+vs1*vp2*rho1*rho2*p4+qa*qa*thetasq*p2*p3*p4)/d;
        
        r[,2 ] = -2.0*vp1*theta*p1*(qa*p3*p4*y+vp2*vs2*x*z)*fac2/d;
        
        r[,3 ]  = 2.0*vp1*rho1*p1*(vs2*p2*x+vs1*p4*y)*fac3/d;

        r[,4] = -2.0*vp1*rho1*theta*p1*(qa*p2*p3-vs1*vp2*z)*fac4/d;

                                        # %%%%%%              /* c    factor to determine rpa(**) */
        fp[1]  = 1.0;
        fp[2]   = vs1/vp1;
        fp[3]  = vp2/vp1;
        fp[4] = vs2/vp1;
        
                                        # %%%%%%              /* c    factor to determine rea(**)    */     
        fe[, 1 ] = 1.0;
        fe[, 2 ] = p2*vp1/(p1*vs1);
        fe[, 3 ] = rho2*p3*vp1/(p1*vp2*rho1);
        fe[, 4] = rho2*p4*vp1/(p1*vs2*rho1)  ;  
      }  else if(incw==2){
                                        #  case 2

        r[, 1 ]  = -2.0*vs1*theta*p2*(qa*p3*p4*y+vp2*vs2*x*z)*fac1/d;
        r[, 2 ]    = 1.0-2.0*p2*(vp2*vs2*p1*x*x+vp1*vs2*rho1*rho2*p3+qa*qa*thetasq*p1*p3*p4)/d;
        
        r[, 3 ]    = 2.0*vs1*rho1*theta*p2*(qa*p1*p4-vp1*vs2*z)*fac3/d;

        r[, 4]   = 2.0*vs1*rho1*p2*(vp1*p3*y + vp2*p1*x)*fac4/d;
        

                                        # %%%%%%              /* c       factor to determine rpa(**) */

        fp[1]  = vp1/vs1;
        fp[2] = 1.0;
        fp[3]   = vp2/vs1 ;     
        fp[4]  = vs2/vs1;

                                        # %%%%%%              /* c   factor to determine rea(**) */
        
        fe[, 1 ] = vs1*p1/(vp1*p2);
        fe[, 2 ] = 1.0;
        fe[, 3 ] = rho2*vs1*p3/(rho1*vp2*p2);
        fe[, 4 ] = rho2*vs1*p4/(rho1*vs2*p2);

        
        
      }


                                        # %%%%%%%%%%%%%%%%%%%%%%%%%
    for(j in 1:4)
      {

        if(icoef==1)
          {
            
                                        # case 1
            rmat[, j ] = Mod(r[, j ]);
          } else if(icoef==2){
                                        # case 2
            rmat[, j] =  Mod(r[, j ]*fp[j]);
          } else if(icoef==3){
                                        # case 3
            rmat[, j] = Mod((r[, j ]*fp[j]) * (r[, j ]*fp[j])  *fe[, j ]);
          }
        

        rra[, j ] = Re(r[, j ]);
        ria[, j ] =  Im(r[, j ]);

        aflag = Mod(rra[, j ])>0.0000001;
        ang[aflag, j ]= atan2(ria[aflag,  j ],rra[aflag, j ]);
        ang[!aflag , j ]= 0.0;
        
        
      }
    return(list(angle=angle, rmat=rmat, rra=rra, ria=ria, ang=ang, incw=incw, icoef=icoef))
  }

