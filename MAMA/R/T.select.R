#Vyber hranicnej hodnoty pre T-statistiku
T.select<-function(stat,fig=TRUE)
{
 quan <- quantile(abs(stat),seq(0.97,0.98,.0001))
 x <- seq(0.97,0.98,1e-04)
if (fig) {
 plot(x,quan,type="b",xlab="percent",ylab="t")
 z <- lm(quan ~ x)
 abline(z,col="red")
 points(0.9787,quan["97.87%"],pch=19,cex=1.5,col="red")
 }
return(quan)
}