corrPlot <-
function(x, y, data,...){
	xx=data[,x]
	yy=data[,y]
	plot(yy~xx,bty="n",mgp=c(2.3,1,0),xlab=x,ylab=y,...)
	lo=loess(yy~xx)
	xs=seq(min(xx),max(xx),l=100)
	lines(xs,predict(lo,newdata=data.frame(xx=xs)),col="red")
	cor=cor.test(xx,yy) 
	mtext(side=3,paste("r =",round(cor$estimate,2),", p =",signif(cor$p.value,2)))
}
