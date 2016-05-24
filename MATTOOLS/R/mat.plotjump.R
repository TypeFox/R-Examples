 
mat.plotjump<-function (inJumpObj) 
{
windows(width=5,height=5)
x=inJumpObj$alphacor[,1]
y=inJumpObj$alphacor[,2]
plot(x,y,xlab="Alpha",ylab="Pearson's R (correlation)",type="l")
par(usr=c(0,1,0,1))
text(0.8,0.8,paste("alpha =",inJumpObj$alpha))

  }



