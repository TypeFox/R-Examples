eqp.eq<-function(sx,sy,X,Ky=max(sy)) 
UseMethod("eqp.eq")

eqp.eq.default<-function(sx,sy,X,Ky=max(sy))
{
	###########################
	#Call parameters
	###########################
	cl<-match.call()

 sx<-sort(sx)
 sy<-sort(sy)
 if (length(X)==1){
   resu<-P.inv(P.x(X,sx),sy,Ky)}
 else{
   temp<-sapply(X,P.x,S=sx)
   resu<-sapply(temp,P.inv,S=sy)
	}
res<-list(call=cl,X=X,resu=resu)
class(res)<-"eqp.eq"
return(res)
  }

print.eqp.eq<-function(x,...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nEquated values under equipercentile Equating:\n")	
	cat("\n")
	print(data.frame(Score=x$X,eqYx=x$resu))
	cat("\n")
}	