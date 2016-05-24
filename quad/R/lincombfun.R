lincombfun <-
function(S,mycoef){

coef.1=mycoef[[1]]
coef.2=mycoef[[2]]
coef.3=mycoef[[3]]
coef.4=mycoef[[4]]

f.1=mycoef[[5]]
f.2=mycoef[[6]]
f.3=mycoef[[7]]
f.4=mycoef[[8]]

 
 # S is a list of sums output from the sum function

 ###########
 ### r=1 ###
 ###########
 r1.p1<-coef.1*S$r1.s1

 ###########
 ### r=2 ###
 ###########
 tmp<-c(S$r2.s1,S$r2.s2) %*% coef.2
 r2.p1<-tmp[1]
 r2.p2<-tmp[2]
 r2.p3<-S$r2.s3-as.numeric(tmp %*% f.2[-length(f.2)])
 
 ###########
 ### r=3 ###
 ###########
 tmp<-c(S$r3.s1,S$r3.s2,S$r3.s3,S$r3.s4,S$r3.s5,S$r3.s6,S$r3.s7) %*% coef.3
 r3.p1<-tmp[1]
 r3.p2<-tmp[2]
 r3.p3<-tmp[3]
 r3.p4<-tmp[4]
 r3.p5<-tmp[5]
 r3.p6<-tmp[6]
 r3.p7<-tmp[7]
 r3.p8<-S$r3.s8-as.numeric(tmp %*% f.3[-length(f.3)])

 ###########
 ### r=4 ###
 ###########
 tmp<-unlist(S[grep("r4",names(S))])[-length(f.4)] %*% coef.4
 r4<-list()
 for(i in 1:(length(f.4)-1)){
  r4[[paste("r4.p",i,sep="")]]<-tmp[i]
  }
 r4[["r4.p23"]]<-S$r4.s23-as.numeric(tmp %*% f.4[-length(f.4)])
   
 #return(list(r1.s1,r2.s1,r2.s2,r2.s3,r3.s1,r3.s2,r3.s3,r3.s4,r3.s5,r3.s6,r3.s7,r3.s8))
 #return(list(r1.p1=r1.p1,r2.p1=r2.p1,r2.p2=r2.p2,r2.p3=r2.p3,r3.p1=r3.p1,r3.p2=r3.p2,r3.p3=r3.p3,r3.p4=r3.p4,r3.p5=r3.p5,r3.p6=r3.p6,r3.p7=r3.p7,r3.p8=r3.p8))
 tmp<-list(r1.p1=r1.p1,r2.p1=r2.p1,r2.p2=r2.p2,r2.p3=r2.p3,r3.p1=r3.p1,r3.p2=r3.p2,r3.p3=r3.p3,r3.p4=r3.p4,r3.p5=r3.p5,r3.p6=r3.p6,r3.p7=r3.p7,r3.p8=r3.p8)
 return(c(tmp,r4))
 }
