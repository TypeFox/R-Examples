YPmodel.setParameter <-
function(startPoint=c(0,0), nm=log(100), maxIter1=50, maxIter2=20, ...)
{
#	if(! is.null(rTestData)){
#		repNum <- rTestData$repNum
#	}

#	if(is.null(Data)){
#		if(is.null(jh)){
#			jh <- 0.001
#		}
#		if(is.null(h)){
#			h <- 0.001
#		}
#	}
#
#	if(! is.null(Data)){
#		if(is.null(jh)){
#			jh <- 1/sqrt(Data$length)
#		}
#		if(is.null(h)){
#			h <- 1/sqrt(Data$length)
#		}
#	}


	Parameters <- list(GroupNum=2, 
		startPoint=startPoint,
		nm=nm, 
		maxIter1=maxIter1,
		maxIter2=maxIter2)

	return(Parameters)
}
