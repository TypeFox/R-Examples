regbelow <-
function(behavior,phaseX, v1,v2){
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  meanA=mean(A,na.rm=T)
  x1=(c(seq(1:tmaxA)))
  
  regA<-lm(A~x1)
  rA<-residuals(regA)
  yA<-regA$coefficients[1]
  BetaA<-regA$coefficient[2]
  
  
  
  tmaxB<-t1[names(t1)==v2]
  startB<-match(v2,phaseX)
  endB<-tmaxB+startB-1
  #tsxB<-behavior[startB:endB]
  B=(behavior[startB:endB])
  x2=(c(seq(1:tmaxB)))
  
  
  cdcl<-c(A,NA,B)
  y<-na.omit(cdcl)
  total=length(y)
  iv=(1:total)
  end<-which(is.na(cdcl))
  iv<-insert(iv,NA,end)
  x2=iv[end+1:total]
  regc<-lm(cdcl~iv)
  
  
  
  x2<-na.omit(x2)
  Byhat<-yA+iv*BetaA
  Byhat<-na.omit(Byhat)
  yhatA<-Byhat[startA:endA]
  startB<-startB-1
  endB<-endB-1
  yhatB<-Byhat[startB:endB]
  maxy=which.max(cdcl)
  
  
  dzone<- y< Byhat
  len1=length(A)
  len2=length(B)
  
  pA<-rep(v1,len1)
  pB<-rep(v2,len2)
  
  p<- c(pA,pB)
  
 
  tm<-table(dzone,p)
  

  
  ctbl<-cbind(tm[,v1],tm[,v2])
  print(ctbl)
  print(prop.table(ctbl,1)*100)
  print(prop.table(ctbl,2)*100)
  c1<-chisq.test(ctbl,correct=FALSE)
 f1<-fisher.test(ctbl,alternative = "two.sided")
  print(c1)
  print(f1)
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  meanA=mean(A,na.rm=T)
  x1=(c(seq(1:tmaxA)))
  
  regA<-lm(A~x1)
  rA<-residuals(regA)
  yA<-regA$coefficients[1]
  BetaA<-regA$coefficient[2]
  
  
  
  tmaxB<-t1[names(t1)==v2]
  startB<-match(v2,phaseX)
  endB<-tmaxB+startB-1
  #tsxB<-behavior[startB:endB]
  B=(behavior[startB:endB])
  
  cdcl<-c(A,NA,B)
  y<-na.omit(cdcl)
  total=length(y)
  iv=(1:total)
  end<-which(is.na(cdcl))
  iv<-insert(iv,NA,end)
  x2=iv[end+1:total]
  regc<-lm(cdcl~iv)
  
  
  
  x2<-na.omit(x2)
  Byhat<-yA+x2*BetaA
  Byhat<-na.omit(Byhat)
  yhatA<-Byhat[startA:endA]
  startB<-startB-1
  endB<-endB-1
  yhatB<-Byhat[startB:endB]
  maxy=which.max(cdcl)
  
  max<-cdcl[maxy]+1
  numx<-sum(!is.na(cdcl))+3
  par(mfrow=c(3,3)) 
  
  maxy=which.max(behavior)
  
  max<-behavior[maxy]+1
  
  numx<-sum(!is.na(behavior))+3
  
  
  graphics.off()
  
  layout(rbind(1,2), heights=c(6,1))
  plot(iv,cdcl, ylim=c(0,max),lwd=2,type="o",col="red", bty="l",xlab="time", ylab="behavior", main="Regression Line" )
  
  abline(reg=regA,col='Blue',lty="dashed")
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("regression line"), lty=c("dashed"), col = c("blue"), lwd = 1,ncol=2,bty ="n") 
  
}
