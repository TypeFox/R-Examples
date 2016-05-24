depth<-function(U,V,m){
#
#  Compute the halfspace depth of the point (u,v) for the pairs of points
#  in the n by 2 matrix m.
#
X<-m[,1]
Y<-m[,2]
FV<-NA
NUMS<-0
NUMH<-0
SDEP<-0.0
HDEP<-0.0
N<-length(X)
P<-acos(-1)
P2<-P*2.0
EPS<-0.000001
ALPHA<-NA
NT<-0
for(i in 1:nrow(m)){
               DV<-sqrt(((X[i]-U)*(X[i]-U)+(Y[i]-V)*(Y[i]-V)))
              if (DV <= EPS){
              NT<-NT+1
                              }
          else{
              XU<-(X[i]-U)/DV
              YU<-(Y[i]-V)/DV
              if (abs(XU) > abs(YU)){
                  if (X[i] >= U){
                      ALPHA[i-NT]<-asin(YU)
                      if(ALPHA[i-NT] < 0.0)
                          ALPHA[i-NT]<-P2+ALPHA[i-NT]
                                  }
                  else{
                      ALPHA[i-NT]<-P-asin(YU)
                      }
                                    }
              else{
                  if (Y[i] >= V)
                      ALPHA[i-NT]<-acos(XU)
                  else
                      ALPHA[i-NT]<-P2-acos(XU)
                  }
              if (ALPHA[i-NT] >= P2-EPS) ALPHA[i-NT]<-0.0
                }
}
NN<-N-NT
if(NN<=1){
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      return(HDEP)
}
ALPHA<-sort(ALPHA[1:NN])
ANGLE<-ALPHA[1]-ALPHA[NN]+P2
for(i in 2:NN){
ANGLE<-max(c(ANGLE,ALPHA[i]-ALPHA[i-1]))
               }
if(ANGLE > (P+EPS)){
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      return(HDEP)
                  }
ANGLE<-ALPHA[1]
NU<-0
for (i in 1:NN){
ALPHA[i]<-ALPHA[i]-ANGLE
if(ALPHA[i]<(P-EPS))NU<-NU+1
               }
if(NU >= NN){
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      return(HDEP)
}
#
#  Mergesort the alpha with their antipodal angles beta,
#  and at the same time update I, F(I), and NBAD.
#
JA<-1
JB<-1
      ALPHK<-ALPHA[1]
      BETAK<-ALPHA[NU+1]-P
      NN2<-NN*2
      NBAD<-0
      I<-NU
      NF<-NN
for(J in 1:NN2){
           ADD<-ALPHK+EPS
          if (ADD < BETAK){
              NF<-NF+1
              if(JA < NN){
                  JA<-JA+1
                  ALPHK<-ALPHA[JA]
              }
              else
                  ALPHK<-P2+1.0
              }
          else{
              I<-I+1
              NN1<-NN+1
              if(I==NN1){
                  I<-1
                  NF<-NF-NN
              }
              FV[I]<-NF
              NFI<-NF-I
              NBAD<-NBAD+depths1(NFI,2)
              if(JB < NN){
                  JB<-JB+1
                  if(JB+NU <= NN)
                      BETAK<-ALPHA[JB+NU]-P
                  else
                      BETAK<-ALPHA[JB+NU-NN]+P
              }
              else
                  BETAK<-P2+1.0
          }
}
NUMS<-depths1(NN,3)-NBAD
#
#  Computation of NUMH for halfspace depth.
#
      GI<-0
      JA<-1
      ANGLE<-ALPHA[1]
      dif<-NN-FV[1]
      NUMH<-min(FV[1],dif)
for(I in 2:NN){
          AEPS<-ANGLE+EPS
          if(ALPHA[I] <= AEPS){
              JA<-JA+1
                              }
          else{
              GI<-GI+JA
              JA<-1
              ANGLE<-ALPHA[I]
              }
          KI<-FV[I]-GI
          NNKI<-NN-KI
          NUMH<-min(c(NUMH,min(c(KI,NNKI))))
   }
NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
depths1(NT,3)
      if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
      NUMH<-NUMH+NT
      HDEP<-(NUMH+0.0)/(N+0.0)
      HDEP
}
