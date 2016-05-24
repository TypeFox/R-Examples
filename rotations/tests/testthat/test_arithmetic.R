rs<-rvmises(20)
Rs<-genR(rs,space='SO3')
Qs<-as.Q4(Rs)
RsRs<-Rs+Rs

RsRs2<-matrix(NA,20,9)
Rst<-matrix(NA,20,9)
Qst<-matrix(NA,20,4)

rs2<-2*abs(rs)

for(i in 1:20){
  
  Rsi<-matrix(Rs[i,],3,3)
  RsRs2[i,]<-Rsi%*%Rsi
  Rst[i,]<-t(Rsi)
  Qst[i,]<-matrix(Qs[i,],nrow=1)*c(1,-1,-1,-1)
  
  if(rs2[i]>pi){
    rs2[i]<-2*pi-rs2[i]
  }
  
}
RsRs2<-as.SO3(RsRs2)
Rst<-as.SO3(Rst)
Qst<-as.Q4(Qst)
context("Arithmetic")

expect_equal(rs2,mis.angle(Rs+Rs))
expect_equal(RsRs,RsRs2)
expect_equal(Rs-Rs,as.SO3(matrix(diag(3),20,9,byrow=T)))
expect_true(all.equal(mean(Rs-mean(Rs)),id.SO3,check.attributes=FALSE))
expect_equal(Rs+Rs-Rs,Rs)
expect_equal(-Rs,Rst)

expect_equal(Qs+Qs,as.Q4(Rs+Rs))
expect_equal(Qs-Qs,as.Q4(matrix(c(1,0,0,0),20,4,byrow=T)))
expect_equal(Qs+Qs-Qs,Qs)
expect_equal(-Qs,Qst)
