YPmodel.estimate <-
function(data=c(), startPoint=c(0,0), nm=log(100), maxIter1=50, maxIter2=20, interval=1, Internal=c(), ...)
{

	if(is.null(data)){
		stop(paste(fun.errorMessage('DataSet')))
	}

	Data <- YPmodel.inputData(data=data)

	Parameters <- YPmodel.setParameter(startPoint=startPoint, nm=nm, maxIter1=maxIter1, maxIter2=maxIter2)


#	if(is.null(Parameters)){
#		warning(paste(fun.errorMessage('DefaultParameter')))
#		Parameters <- YPmodel.setParameter(Data=Data)
#	}

#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#

	maxIteration1 <- Parameters$maxIter1
	maxIteration2 <- Parameters$maxIter2
	nm <- Parameters$nm
	startPoint <- Parameters$startPoint

	Z <- Data$Z
	Delta <- Data$Delta
	n <- Data$length

	nb <- array(startPoint,c(1,2))
	ob <- nb + 0.1
	bc <- t(nb)

	data2 <- fun.hcvxitr1(bc,Data,Parameters)
	po <- data2$po
	s <- data2$s

	of <- s[1]
	#nm <- 100

	numIteration <- 0

	#-----------------------------------------------------------------#
	while ((sum(abs(round(1000*ob) - round(1000*nb)))>0) & (numIteration < maxIteration1) ) {

	    numIteration <- numIteration + 1
	    t <- 1
	    bc <- t(nb) + t*po
	    idb <- sum((bc > nm) + (bc < -nm))

		numIterationInner<- 0

	   #-----------------------------------------------------------------#
	   while ((idb > 0) & (numIterationInner < maxIteration2)){
	        numIterationInner<- numIterationInner + 1
	        t <- t*0.5
	        bc <- t(nb)+ t*po
	        idb <- sum(( bc > nm) + (bc < - nm))
	    }
	   #-----------------------------------------------------------------#

	    m1 <- 1
	    b <- bc

	   ##return [s,ru]
	   #data4 <- fun.ntitr(b,Z,Delta,n,m1)
	   #s <- data4$s
	   #ru <- data4$ru

		data4 <- fun.oldp2(b,1,Data)
		s <- data4$s
		ru <- data4$ru

	    if (of==0){
	        break
	    } else if (abs((s/of-1)*1.e+8)<1){
	        break
	    }

	    ob <- nb
	    nb <- t(bc)
	    of <- s

	   ##return [po,s]
	   #data3 <- fun.hcvxitr(bc,jh,Z,Delta,n)
	   #po <- data3$po
	   #s <- data3$s

	   data3 <- fun.hcvxitr1(bc, Data, Parameters)
	   po <- data3$po
	   s <- data3$s

	}
	#-----------------------------------------------------------------#

	b <- cbind(t(nb), t(ob))
	m <- 2

	##return [s,ru]
	#data5 <- fun.ntitr(b,Z,Delta,n,m)
	#s <- data5$s
	#ru <- data5$ru

	data5 <- fun.oldp2(b,m,Data)
	s <- data5$s
	ru <- data5$ru

	ots <- min(s)
	jt <- which.min(s)
	beta <- t(b[,jt])
	r <- ru[,jt]

	Estimate<- list(beta=beta,r=r)

	#-----------------------------------------------------------------#
	if(interval==1){

		if(is.null(Internal)){
			Internal <- fun.internalParameters(Data, Estimate)
		}

		ru <- Internal$ru
		p <- Internal$p
		pl <- Internal$pl
		bt <- Internal$bt
		deni <- Internal$deni
		kall <- Internal$kall
		sm <- Internal$sm
		gama <- Internal$gama
		b <- Internal$b

		pq <-fun.variance(b,bt,ru,gama,p,pl,deni,sm,kall,Data)
		variance.beta1<-pq[1,1]
		variance.beta2<-pq[1,2]

	#-----------------------------------------------------------------#
	## Output Resuts
	#-----------------------------------------------------------------#
	    Estimate$variance.beta1 <- variance.beta1
	    Estimate$variance.beta2 <- variance.beta2

	    Estimate$Data <-Data
	    Estimate$Parameters <- Parameters
	    Estimate$interval <- interval

	#-----------------------------------------------------------------#
	}

	class(Estimate) <- "YPmodel.estimate"
	Estimate$call <- match.call()

	return(Estimate)

}
