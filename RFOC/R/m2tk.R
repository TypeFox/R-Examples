m2tk<-function(m0)
{
   #m consists of diagonal components of a diagonalized moment tensor
   ###  Author Keehoon Kim, UNC, 2012
   mo=sort(m0,decreasing=TRUE)
   M=sum(mo)/3
   Md=mo-M
   if(Md[2]>=0) {k=M/(abs(M)-Md[3])}
   if(Md[2]<0) {k=M/(abs(M)+Md[1])}
   if(Md[2]>0) {T=-2*Md[2]/Md[3]} 
   if(Md[2]==0) {T=0}
   if(Md[2]<0) {T=2*Md[2]/Md[1]}
   return(list(k=k,T=T))
}
