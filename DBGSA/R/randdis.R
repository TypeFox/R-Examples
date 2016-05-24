randdis<-function(z,minigenenum,randnum,setdis,normnum,Meth,resultname){
z<- as.matrix(z)
result<-matrix( ,minigenenum,randnum);

for(n in 1:minigenenum)
 {
  for(m  in 1:randnum)
     {
      
      A<-sample(1:nrow(z),n);#产生随机数

      B<-matrix(,n,ncol(z))
      for(i in 1:n )
         {
         
         B[i,]=z[A[i],];
         }
      
         result[n,m]=setdis(B,normnum,Meth)
           }
   }
write.table(result,file=resultname,quote=FALSE)
}




