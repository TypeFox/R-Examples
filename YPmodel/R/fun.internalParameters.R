fun.internalParameters <-
function(Data, Estimate=c())
{

	beta <- Estimate$beta
	r <- Estimate$r

	#return [fb1,fb2]
	data6<-fun.generateParameters(Data)
	# 
	GroupData <- data6$GroupData
	fb <- data6$fb
	#K(t) = \sum_{i=1}^n I(X_i \geq t)#
	kall <- data6$kall
	u0 <- data6$u0
	kk <- data6$kk
	l <- data6$l
	dfb1 <- data6$dfb1

	bt <- exp(-c(beta[1],beta[2]))
	b <- t(beta)

	data7 <- fun.oldp2(b,1,Data)
	s <- data7$s
	#
	ru <- data7$ru
	u <- data7$u
	#
	gama <- data7$gama 
	p <- data7$p
	pl <- data7$pl 
	deni <- data7$deni 
	sm <- data7$sm 

	InternalParameters <- list(
		beta=beta,
		r=r,
	  	ru=ru,
	  	fb=fb,
		bt=bt,
	  	kk=kk,
	  	p=p,
	  	pl=pl,
	  	deni=deni,
	  	kall=kall,
	  	sm=sm,
	  	gama=gama,
	  	b=b,
	  	u0=u0,
	  	l=l,
	  	dfb1=dfb1,
	  	GroupData=GroupData)

	return(InternalParameters)

}
