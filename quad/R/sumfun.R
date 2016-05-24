sumfun <-
function(W){

 # useful
 W2<-W^2
 rowsums.W<-rowSums(W)
 rowsums.W2<-rowSums(W^2)
 rowsums.W3<-rowSums(W^3)
 rowsums.W4<-rowSums(W^4)
 oW<-outer(rowsums.W,rowsums.W)
 sum.W<-sum(rowsums.W)
 sum.W2<-sum(rowsums.W2)
 sum.W3<-sum(rowsums.W3)
 sum.W4<-sum(rowsums.W4)
 WW<-W %*% W
 
 ###########
 ### r=1 ###
 ###########
 r1.s1<-sum.W

 ###########
 ### r=2 ###
 ###########
 r2.s1<-sum.W2
 r2.s2<-sum(rowsums.W^2)	
 r2.s3<-r1.s1^2

 ###########
 ### r=3 ###
 ###########
 r3.s1<-sum.W3
 r3.s2<-sum(rowsums.W*rowsums.W2)
 r3.s3<-sum(W*oW)
 r3.s4<-sum(W*WW)					
 r3.s5<-sum(rowsums.W^3)
 r3.s6<-r1.s1*r2.s1
 r3.s7<-r1.s1*r2.s2
 r3.s8<-r1.s1^3

 ###########
 ### r=4 ###
 ###########

 r4<-list()
 r4[["r4.s1"]]<-sum.W4
 r4[["r4.s2"]]<-sum(rowsums.W3*rowsums.W)
 r4[["r4.s3"]]<-sum(rowsums.W2^2)
 r4[["r4.s4"]]<-sum(W2*WW)
 r4[["r4.s5"]]<-sum(rowsums.W2*rowsums.W^2)
 r4[["r4.s6"]]<-sum(W2*oW)
 r4[["r4.s7"]]<-sum(W*outer(rowsums.W2,rowsums.W))
 r4[["r4.s8"]]<-sum(WW^2)
 r4[["r4.s9"]]<-sum(rowsums.W*rowSums(W*WW))
 r4[["r4.s10"]]<-sum((W %*% rowsums.W)^2)
 r4[["r4.s11"]]<-sum(rowsums.W^4)
 r4[["r4.s12"]]<-sum(W*outer(rowsums.W^2,rowsums.W))
 r4[["r4.s13"]]<-r2.s1^2
 r4[["r4.s14"]]<-r2.s2^2
 r4[["r4.s15"]]<-r2.s1*r2.s2
 r4[["r4.s16"]]<-r2.s1*r2.s3
 r4[["r4.s17"]]<-r2.s2*r2.s3
 r4[["r4.s18"]]<-r1.s1*r3.s1
 r4[["r4.s19"]]<-r1.s1*r3.s2
 r4[["r4.s20"]]<-r1.s1*r3.s3
 r4[["r4.s21"]]<-r1.s1*r3.s4
 r4[["r4.s22"]]<-r1.s1*r3.s5
 r4[["r4.s23"]]<-r1.s1^4

 #return(list(r1.s1,r2.s1,r2.s2,r2.s3,r3.s1,r3.s2,r3.s3,r3.s4,r3.s5,r3.s6,r3.s7,r3.s8))
 #return(list(r1.s1=r1.s1,r2.s1=r2.s1,r2.s2=r2.s2,r2.s3=r2.s3,r3.s1=r3.s1,r3.s2=r3.s2,r3.s3=r3.s3,r3.s4=r3.s4,r3.s5=r3.s5,r3.s6=r3.s6,r3.s7=r3.s7,r3.s8=r3.s8))
 tmp<-list(r1.s1=r1.s1,r2.s1=r2.s1,r2.s2=r2.s2,r2.s3=r2.s3,r3.s1=r3.s1,r3.s2=r3.s2,r3.s3=r3.s3,r3.s4=r3.s4,r3.s5=r3.s5,r3.s6=r3.s6,r3.s7=r3.s7,r3.s8=r3.s8)
 return(c(tmp,r4))
 }
