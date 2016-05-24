cubic=function(A,B,C,D){
### Function that gets the 4 variables of the cubic equation startuing with the cubic x factor as A.
### All the values are real number.
e1=complex(1);
e2=complex(1);
a=numeric(1);
b=numeric(1);
c=numeric(1);
d=numeric(1);
p=numeric(1);
q=numeric(1);
Dis=complex(1);
u=complex(1);
v=complex(1);
gou=numeric(1);
gov=numeric(1);
test=complex(1);
z=c();
### It is checked wheter only a contstant is remaining.
### In this case the an empty numeric vector is returned.
if(A==0 & B==0 & C==0){
  z[1]=0;
  z[2]=0;
  z[3]=0;
  warning("ABC all zero should not be possible when e is not zero\n");
  return(z);
  };
### If the firts two factors are zero only a linear function remains.
### The solution is simple, because C is not zero. This case has been selected in the 
### condition before.
if(A==0 & B==0){
   z[1]=(-D)/C;
   z[2]=z[1];
   z[3]=z[2];
   return(z);
   };
### If the equatio is of degree 2 a closed solution is known. 
### Please recognoze that B can not be zero becuase this case is solved above already.
if(A==0){
   z[1]=(-C/(2*B))+sqrt((C/(2*B))**2-(D/B));
   z[2]=(-C/(2*B))-sqrt((C/(2*B))**2-(D/B));
   z[3]=z[2];
   return(z);
   }
### The coefficients are normalized to get a 1 infront of the cubic x.
### The coefficients are renamed to fit the names of the solution formular.
### Now a is infront of x square and d is not used any more.
d=A/A;
a=B/A;
b=C/A;
c=D/A;
### Linear transformation to get rid of the squared coefficient.
p=b-(1/3)*(a**2);
q=(2/27)*(a**3)-(1/3)*(a*b)+c;
Dis=(q/2)**2+(p/3)**3;
e1=-0.5+0.5*1i*sqrt(3);
e2=-0.5-0.5*1i*sqrt(3);
u=(-0.5*q+(sqrt(Dis+0i)))**(1/3);
v=(-0.5*q-(sqrt(Dis+0i)))**(1/3);
#print("Startwerte u und v und Dis");
#print(u);
#print(v);
#print(Dis);
#print(p);
#print(q);
### All three solutions are printed as a solution multiplied by a cubic unique root.
### Since 3*u*v+p=0 it is u*v=-(1/3)*p which is checked to be close to zero. This 
### necesary for numerical stability.
### All alternativies are tested and the variables u and v are 
### turned in the two dimensional space by a cubic untit root 
### until this condition is fullfilled.
###
for(gou in c(1:3)){
  for(gov in c(1:3)){
      test=(u*v+(1/3)*p);
#      print(p);
#      print("Angfang")
#      print(u);
#      print(v);
#      print(test);
#      print("next");
      if(abs(u*v-(-1)*(1/3)*p)<0.00001) break;
      v=v*e1;
    };
  if(abs(u*v-(-1)*(1/3)*p)<0.00001) break;
  u=u*e1;
};
#### t is the solution. Stored in z elements one two and three.
####
z[1]=u+v;
z[2]=u*e1+v*e2;
z[3]=u*e2+v*e1;
### Now the linear transformation is backwards used. Since a is b and the coefficent of cubic x is 1 t_untransposed=t -(a/3*d) and d=1. 
### This leads numeric stable to the cubic roots.
### The solutions are eighter of rank one or rank two and complex roots are symmetric. 
###
### linear Backtranformation after normalization of the cubic x coefficient
z=z-a/3;
z[10]=u*e1+v*e1;
z[11]=u*e2+v*e2;
z[4]=p;
z[5]=q;
z[6]=Dis;
return(z);
}
optimalfrequency=function(mm,mmorg){
### The mm is the 2x2 matrix and mmorg is the corresponding 3x3 matrix of the unsorted haplotypes.
### The result is an object of the class "optimalfrequency". The attributes are T value of the t distribution named "Testvalue".
### The propability of not being zero by change called: "prSimilarByChange"
D=det(mm);
r1=sum(mm[,1]);
r2=sum(mm[,2]);
r3=sum(mm[1,]);
r4=sum(mm[2,]);
n=sum(as.numeric(mm));
n2=sum(as.numeric(mmorg));
mmorg=matrix(as.numeric(as.matrix(mmorg))*(1/n2),nrow=length(mmorg[,1]));
LK=D/sqrt(r1*r2*r3*r4);
r11=r3/n;
r12=r4/n;
s11=r1/n;
s12=r2/n;
#n2=sum(as.numeric(mmorg));
rest=(((r11*r11*s11*s11-as.numeric(mmorg[1,1]))**2)/(r11*r11*s11*s11))+(((r11*r11*s11*s12-as.numeric(mmorg[1,2]))**2)/(r11*r11*s11*s12))+
     (((r11*r11*s12*s12-as.numeric(mmorg[1,3]))**2)/(r11*r11*s12*s12))+
     (((r11*r12*s11*s11-as.numeric(mmorg[2,1]))**2)/(r11*r12*s11*s11))+(((r11*r12*s11*s12-as.numeric(mmorg[2,2]))**2)/(r11*r12*s11*s12))+
     (((r11*r12*s12*s12-as.numeric(mmorg[2,3]))**2)/(r11*r12*s12*s12))+
     (((r12*r12*s11*s11-as.numeric(mmorg[3,1]))**2)/(r12*r12*s11*s11))+(((r12*r12*s11*s12-as.numeric(mmorg[3,2]))**2)/(r12*r12*s11*s12))+
     (((r12*r12*s12*s12-as.numeric(mmorg[3,3]))**2)/(r12*r12*s12*s12));
rest2=rest*(n2); # theoretisch falsch tested chiquadrat verteilung gibt p wert von 1
#rest2=rest;
pr=pchisq(rest2,8);
pr=1-pr;
result=list();
class(result)="optimalfrequency";
result[["LK"]]=LK;
if(r1*r2*r3*r4==0){rest2=-999999999999999999;}
result[["Testvalue"]]=rest2;
result[["prSimilarByChange"]]=pr;
#return(LK,rest2,pr);
return (result);
}
findoptimal=function(A,B,C,D,mmorg,exact=0.00001){
### This programs gets the four element of the solution haplotype matrix and the original 3x3 matrix mmorg
###
###
###
goall=c();
ppall=list();
rest=c();
LK=list();
FitGo=c();
### limiting the number of results that will be returned or coerce to to a numeric contstant
genau=numeric(1);
AA=numeric(3);
BB=numeric(3);
CC=numeric(3);
DD=numeric(3);
genau=exact;
AA=-1;
BB=-1;
CC=-1;
DD=-1;
mm=matrix(rep(0,4),nrow=2)
pp=numeric(1);
p=numeric(3);
p=c();
pmax=numeric(1);
zeahler=c();
maxzaehler=numeric(1);
#class(result)=c("2 Snp with two cubic predicte haplotype counts");
for (go in c(1:length(A))){
### limit is exact paramter
### Handling the real part as almost zero
if(abs(Re(A[go]))<genau){A[go]=A[go]-Re(A[go]);};
if(abs(Re(B[go]))<genau){B[go]=B[go]-Re(B[go]);};
if(abs(Re(C[go]))<genau){C[go]=C[go]-Re(C[go]);};
if(abs(Re(D[go]))<genau){D[go]=D[go]-Re(D[go]);};   
#  
### Handling the imaginary part as almost zero 
if(abs(Im(A[go]))<genau) { A[go]=Re(A[go])};
if(abs(Im(B[go]))<genau) { B[go]=Re(B[go])};
if(abs(Im(C[go]))<genau) { C[go]=Re(C[go])};
if(abs(Im(D[go]))<genau) { D[go]=Re(D[go])};
if(Im(A[go])==0){ AA[go]=A[go]};
if(Im(B[go])==0){ BB[go]=B[go]};
if(Im(C[go])==0){ CC[go]=C[go]};
if(Im(D[go])==0){ DD[go]=D[go]};
# if all values are negative by zero within on solution than the variables are set to zero
#
### If no sensfully result is avaible this case will be left. For example the cubic polynom contains only D.
if(Re(AA[go])>=0 & Re(BB[go])>=0 &Re(CC[go])>=0 & Re(DD[go])>=0 & ! is.na(AA[go]+BB[go]+CC[go]+DD[go])){
  mm=matrix(c(Re(AA[go]),Re(BB[go]),Re(CC[go]),Re(DD[go])),nrow=2);
#  ppall[[go]]=chisq.test(mm,simulate.p.value=TRUE,B=1000);
#  pp=chisq.test(mm,simulate.p.value=TRUE,B=1000)$p.value;
#  pp=ppall[[go]]$p.value;
  ppall[[go]]=1;
  pp=1;
  p=c(p,pp);
  goall=c(goall,go)
  pmax=max(p);
  zaehler=go;
#  maxzaehler=goall[which(p==pmax)];
#  print(AA);
  LL=list();
  class(LL)="optimalfrequency";
  LL=optimalfrequency(mm,mmorg);
  class(LK[go])="optimalfrequency";
  LK[[go]]=LL;
  FitGo[go]=LK[[go]]$Testvalue;
  maxzaehler=max(maxzaehler,which(abs(LK[[go]]$Testvalue)==min(abs(FitGo),na.rm=TRUE))*go,na.rm=TRUE);
#  rest[go]=LL[[2]];
  };
};
#return (pmax,p,A,B,C,D,zaehler,maxzaehler);
output=list(AA[maxzaehler],BB[maxzaehler],CC[maxzaehler],DD[maxzaehler],pmax,p,LK,ppall,goall)
names(output)=c("AA","BB","CC","DD","pmax","p","LK","ppall","goall")
#return(AA[maxzaehler],BB[maxzaehler],CC[maxzaehler],DD[maxzaehler],pmax,p,LK,ppall,goall);
return(output)
}
haplotypeit=function(a,b,c,d,e,f,g,h,i){
### initialization of the variables
### AA,BB,CC,DD are the coefficients of the cubic polynom
### A1,B1,C1,D1 are the coefficients of the haplotype 2x2 matrix
A1=numeric(1);
B1=numeric(1);
C1=numeric(1);
D1=numeric(1);
A1=2*a+b+d;
B1=2*c+b+f;
C1=2*g+h+d;
D1=2*i+h+f;
AA=numeric(1);
BB=numeric(1);
CC=numeric(1);
DD=numeric(1);
A=complex(3);
B=complex(3);
C=complex(3);
D=complex(3);
AA=2*e*e;
BB=-e*e+A1*e+D1*e-B1*e-C1*e-2*e*e;
CC=(-1)*A1*e-D1*e+A1*D1+B1*C1+B1*e+C1*e+e*e;
DD=(-1)*A1*D1;
### Shotpath to the solution if the entry e in the 3x3 matrix is zero
if(e>0){
   erg=cubic(AA,BB,CC,DD);
 }else{erg=c();erg[1]=0;erg[2]=0;erg[3]=0;};
A[1]=A1+erg[1]*e;
A[2]=A1+erg[2]*e;
A[3]=A1+erg[3]*e;
D[1]=D1+erg[1]*e;
D[2]=D1+erg[2]*e;
D[3]=D1+erg[3]*e;
B[1]=B1+(1-erg[1])*e;
B[2]=B1+(1-erg[2])*e;
B[3]=B1+(1-erg[3])*e;
C[1]=C1+(1-erg[1])*e;
C[2]=C1+(1-erg[2])*e;
C[3]=C1+(1-erg[3])*e;
#### Four vectors are returned. Since the logic of R cahnged when switching version. 
#### The routine returns a list named A,B,C,D of vectors
####
output=list(A,B,C,D)
names(output)=c("A","B","C","D")
#return(A,B,C,D)
return(output)
}
callhaplotype=function(dd){
### this functions gets a 3x3 data.frame 
### temp3daten is the numeric data.frame of the 3x3 snip1, snip2 crossmatch. 
        tempbootobj=list();
        tempbootci =list();
        haplo=list();
        resultvector=c();
        temp2haplo =as.numeric(t(dd));
        temphaplo  =haplotypeit(temp2haplo[1],temp2haplo[2],temp2haplo[3],temp2haplo[4],temp2haplo[5],temp2haplo[6],temp2haplo[7],temp2haplo[8],temp2haplo[9]);
        haplo=findoptimal(temphaplo[[1]],temphaplo[[2]],temphaplo[[3]],temphaplo[[4]],dd);
        resultvector=Re(matrix(c(haplo[[1]],haplo[[2]],haplo[[3]],haplo[[4]]),ncol=2,nrow=2));
        return(resultvector);
}
