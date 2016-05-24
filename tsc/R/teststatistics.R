teststatistics <-
function(x,y){
  delta<-0.1
  length1<-length(x)
  length2<-length(y)
  Lst=list(x,y)
  comb<-c(x,y) 
  sx<-sort(x) # sorting of the generated paired data z1 
  sy<-sort(y)
  scomb<-sort(comb)
  ############Term A############## 
  # mean of each group k=2 
  avg=mean(comb) 
  # length of each group k=2
  lenVec = unlist(lapply(Lst, FUN="length")) 
  sigmaSq = sum((comb - avg)^2)/length(comb) 
  # computer term A and term B for each k 
  termLogA = (length1/2)*log(2*pi*exp(1)*sigmaSq)
  termLogB = (length2/2)*log(2*pi*exp(1)*sigmaSq)
  
  m1=c(1:min(c(round((length1)^(1-delta)),round(length1/2))))
  a<-replicate(length1,m1)           ###store repeated values of a vector "m1"
  rm<-as.vector(t(a))         ###transpose the previous length(m)*n1 matrix and make it to be a vector
  L<-c(1:length1)- rm              ###order from (1-m) to (n1-m)
  LL<-replace(L, L <= 0, 1 )  ###replace value<=0 with 1 when (1-m) <=0
  U<-c(1:length1)+ rm              ###order from (1+m) to (n1+m)
  UU<-replace(U, U > length1, length1)  ###replace value>n1 with n1 when (n1+m)>n1
  
  xL<-sx[LL]  ###obtain x(i-m)
  xU<-sx[UU]  ###obtain x(i+m)
  F<-xU-xL
  F[F==0]<-1/(length1+length2)
  I<-2*rm/(length1*F)   
  ux<-array(I, c(length1,length(m1)))   ### make previous vector a n1 * length(m) matrix
  tstat1<-log(min(apply(ux,2,prod))) 
  m2=c(1:min(c(round((length2)^(1-delta)),round(length2/2))))
  a1<-replicate(length2,m2)           ###store repeated values of a vector "m1"
  rm1<-as.vector(t(a1))         ###transpose the previous length(m)*n1 matrix and make it to be a vector
  L1<-c(1:length2)- rm1              ###order from (1-m) to (n1-m)
  LL1<-replace(L1, L1 <= 0, 1 )  ###replace value<=0 with 1 when (1-m) <=0
  U1<-c(1:length2)+ rm1              ###order from (1+m) to (n1+m)
  UU1<-replace(U1, U1 > length2, length2)  ###replace value>n1 with n1 when (n1+m)>n1
  
  xL1<-sy[LL1]  ###obtain x(i-m)
  xU1<-sy[UU1]  ###obtain x(i+m)
  F<-xU1-xL1
  F[F==0]<-1/(length1+length2)
  I<-2*rm1/(length2*F)   
  ux1<-array(I, c(length2,length(m2)))   ### make previous vector a n1 * length(m) matrix
  tstat2<-log(min(apply(ux1,2,prod))) 
  T<- termLogA+ termLogB+tstat1+tstat2
  return(T)
}
