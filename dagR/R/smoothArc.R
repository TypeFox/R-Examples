smoothArc <-
function(A, B, C, res=20, gap=0.05, p=FALSE)
{ # this has replaced gassoc() as the function to draw associations;

  angBA<-angle(B,A);
  angBC<-angle(B,C);
  beta1<-inAngle(angBA, angBC);

  lenBD<-cos(beta1)*distPoints(B,C); # cos(beta1)=[BD]/[BC];
  D<-anglePoint(B, angBA, lenBD);

  lenCD<-distPoints(C, D);
  lenBC<-distPoints(B, C);

  lenBM<-(lenCD^2+lenBD^2)/(2*lenCD);

  beta<-acos((lenBC/2)/lenBM);

  angBM<-addAngle(angBC, -sign(beta1)*beta);
  M<-anglePoint(B, angBM, lenBM);

  angMC<-angle(M,C);
  angMB<-angle(M,B);
  stepWidth<-inAngle(angMC, angMB);
   stepWidth<-sign(stepWidth)*(abs(stepWidth)-gap/lenBM);
   stepWidth<-stepWidth/res;

  circ1_ang<-seq(from=angMC, by=stepWidth,
                 length.out=res+1);
  points(matrix(anglePoint(M, circ1_ang, lenBM), ncol=2),
         type='l', lty=2);

  angAB<-angle(A,B);
  angAC<-angle(A,C);
  beta1<-inAngle(angAB, angAC);

  lenAD<-cos(beta1)*distPoints(A,C);
  D<-anglePoint(A, angAB, lenAD);

  lenCD<-distPoints(C, D);
  lenAC<-distPoints(A, C);

  lenAM<-(lenCD^2+lenAD^2)/(2*lenCD);

  beta<-acos((lenAC/2)/lenAM);

  angAM<-addAngle(angAC, -sign(beta1)*beta);
  M<-anglePoint(A, angAM, lenAM);

  angMC<-angle(M,C);
  angMA<-angle(M,A);
  stepWidth<-inAngle(angMC, angMA);
   stepWidth<-sign(stepWidth)*(abs(stepWidth)-gap/lenAM);
   stepWidth<-stepWidth/res;

 circ2_ang<-seq(from=angMC, by=stepWidth,
                 length.out=res+1);
  points(matrix(anglePoint(M, circ2_ang, lenAM), ncol=2),
         type='l', lty=2);

  if(p==TRUE) {
    points(C[1], C[2]);
  }
}

