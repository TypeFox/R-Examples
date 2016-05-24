mijsdr<-function(mxx,myy,mzz,mxy,mxz,myz)
  {
    
#############   convert CMT moment tensor form of focal mech to
#############   strike dip rake format
#####   INPUT
#####	mij - siz independent components of the moment tensor
#####
#####   OUTPUT
#####	str - strike of first focal plane (degrees)
#####	dip - dip of first focal plane (degrees)
#####	rake - rake of first focal plane (degrees)
#####
########   extracted from matlab code
########           http://www.ceri.memphis.edu/people/olboyd/Software/bb.m
######### Adapted from matlab code mij2sdr.m and from code, mij2d.f,
########   created by Chen Ji and given to me by Gaven Hayes.  


    ###########  internal function called by main code
    TDL<-function(AN,BN)
      {
##########  AN and BN are vectors

        XN=AN[1];
        YN=AN[2];
        ZN=AN[3];
        XE=BN[1];
        YE=BN[2];
        ZE=BN[3];
        
        AAA=1.0E-06;
        CON=57.2957795;
        
        if (abs(ZN) < AAA){
          FD=90.;
          AXN=abs(XN);
          if (AXN > 1.0) { AXN=1.0;}

          FT=asin(AXN)*CON;
          ST=-XN;
          CT=YN;
          if (ST >= 0. & CT < 0)  { FT=180.-FT; }
          
          if (ST < 0. & CT <= 0)  {   FT=180.+FT; }
          
          if (ST < 0. & CT > 0) {   FT=360.-FT; }
          
          FL=asin(abs(ZE))*CON;
          SL=-ZE;
          if (abs(XN) < AAA) { CL=XE/YN; }
          else
            {
              CL=-YE/XN;
            }

          if (SL >= 0. & CL < 0) {  FL=180.-FL; }

          if (SL < 0. & CL <= 0)  { FL=FL-180.; }

          if (SL < 0. & CL > 0) {   FL=-FL; }
        }
        else
          {
            if (-ZN > 1.0) {  ZN=-1.0; }

            FDH=acos(-ZN);
            FD=FDH*CON;
            SD=sin(FDH);
            if  (SD == 0){ return(list(ft=FT,fd=FD,fl=FL)) }

            ST=-XN/SD;
            CT=YN/SD;
            SX=abs(ST);
            if (SX > 1.0) 
              { SX=1.0; }

            FT=asin(SX)*CON;
            if (ST >= 0. & CT < 0)  { FT=180.-FT; }

            if (ST < 0. & CT <= 0)  {   FT=180.+FT; }

            if (ST < 0. & CT > 0) {   FT=360.-FT; }

            SL=-ZE/SD;
            SX=abs(SL);
            if (SX > 1.0) {    SX=1.0; }
            
            FL=asin(SX)*CON;
            if (ST == 0)
              {  CL=XE/CT; }
            else {
              XXX=YN*ZN*ZE/SD/SD+YE;
              CL=-SD*XXX/XN;
              if (CT == 0) {   CL=YE/ST; }

            }
            if (SL >= 0. & CL < 0) 
              {  FL=180.-FL; }

            if (SL < 0. & CL <= 0) {  FL=FL-180. }

            if (SL < 0. & CL > 0) 
              {  FL=-FL; }

            
          }

        return(list(ft=FT,fd=FD,fl=FL))
        

      }

    


    
    a = matrix(c(mxx, mxy, mxz, mxy, myy, myz, mxz, myz, mzz), ncol=3)
    E = eigen(a)
    
    D = E$values[c(1,3,2) ]
    
    V  = E$vectors
###  for some reason the matlab code does this rearrangement...
    V = -V[, 3:1]

    
    V[2:3,1:3] = -V[2:3,1:3];
    
    V = matrix( c(V[2,3], V[2,1], V[2,2], V[3,3], V[3,1], V[3,2], V[1,3], V[1,1], V[1,2]), ncol=3, byrow=TRUE) ;



    
    IMAX = which(D == max(D));
    IMIN = which(D == min(D));


    AE = (V[,IMAX]+V[,IMIN])/sqrt(2.0);
    AN = (V[,IMAX]-V[,IMIN])/sqrt(2.0);
    
    AER = sqrt(AE[1]^2+AE[2]^2+AE[3]^2);
    ANR = sqrt(AN[1]^2+AN[2]^2+AN[3]^2);
    AE = AE/AER;
    AN = AN/ANR;

    
    if (AN[3] <= 0.)
      {
	AN1 = AN;
	AE1 = AE;
      } else {
	AN1 = -AN;
	AE1 = -AE;
      }



###   [ft,fd,fl]
    FTDL = TDL(AN1,AE1);
    str = 360 - FTDL$ft;
    dip = FTDL$fd;
    rake = 180 - FTDL$fl;

    return(c(str, dip, rake ))
  }






