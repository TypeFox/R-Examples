mea.eq<-function(sx,sy,scale) 
UseMethod("mea.eq")

mea.eq.default<-function(sx,sy,scale)
{
	###########################
	#Call parameters
	###########################
	cl<-match.call()

	mu.x<-mean(sx)
	mu.y<-mean(sy)
	resu<-scale-mu.x+mu.y
	res<-list(call=cl,scale=scale,resu=resu)
	class(res)<-"mea.eq"
	return(res)
}

print.mea.eq<-function(x,...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nEquated values under mean Equating:\n")	
	cat("\n")
	print(data.frame(Score=x$scale,eqYx=x$resu))
	cat("\n")
}	