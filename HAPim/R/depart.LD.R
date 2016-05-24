`depart.LD` <-
function(perf,CD){

 CDcarre	=	CD*CD
 nbre.desc.tot  = 	length(perf)
 mu    		= 	sum(perf*CDcarre)/sum(CDcarre)
 Y     		= 	perf-mu
 sigma 		= 	sum(Y^2*CDcarre)/nbre.desc.tot

 param		=	c(mu,sigma)

 param


}

