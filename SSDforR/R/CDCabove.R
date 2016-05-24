CDCabove <-
function(behavior,phaseX,v1,v2){
  
  nsucces=c(0,0,0,0,5,6,6,7,8,8,9,9,10,11,12,12,12,13,13,14,14,15,15,17,17,18,18,19,19,20,20,21,22,23,23,24,24,25,25,26,27,27,28,28,29,29,30,31,31,32)
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  meanA=mean(A,na.rm=T)
  sdA=(sd(A,na.rm=T))*.25
  meanA<-meanA+sdA
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
  yA<-yA+sdA
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
  miny=which.min(behavior)
  min<-behavior[miny]
  
  numx<-sum(!is.na(behavior))+3
  
  
  
 
  
  
  naboveline<-B>Byhat
  nabovemean<-B>meanA
  
  nabove<-table(naboveline,nabovemean)
  
  needed=nsucces[tmaxB]
  
  lin1<-c("needed=", needed)
  print(lin1)

  writeLines("-----------------above lines------------------")
 writeLines ("TRUE, TRUE = Number above the lines")
  
  print(nabove)
  
  
  
  graphics.off()
  
  layout(rbind(1,2), heights=c(6,1))
  plot(iv,cdcl, ylim=c(min,max),lwd=2,type="o",col="red", bty="l",xlab="time", ylab="behavior", main="CDC" )
  
  abline(h=(meanA),col="green")
  abline(a=yA,b=BetaA,col='Blue',lty="dashed")
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("Adj. regression line","Adj. mean line"), col = c("blue","green"), lwd = 1,ncol=2,bty ="n")  
 
}
