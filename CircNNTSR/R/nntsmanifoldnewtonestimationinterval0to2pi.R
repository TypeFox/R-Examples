nntsmanifoldnewtonestimationinterval0to2pi <-
function(data,cutpoints,subintervals,M=0,iter=1000,initialpoint=FALSE,cinitial)
{
	frequencies <- data/sum(data)
	subintervals<-length(data)

	statisticsarray<-array(0,c(M+1,M+1,subintervals))

	for (k in 1:subintervals){
	     for (k1 in 0:M){
		for (k2 in 0:M){
		     if (k1 != k2)
			statisticsarray[k1+1,k2+1,k] <- -Conj((1i/(k1-k2))*(exp(1i*(k1-k2)*cutpoints[k+1]) - exp(1i*(k1-k2)*cutpoints[k])))
		     else
			statisticsarray[k1+1,k2+1,k] <- (cutpoints[k+1] - cutpoints[k]) 
		}
	     }
	}


	if (initialpoint)
		c0<-cinitial
	else
		c0<-nntsrandominitial(M)

	c0<-c0/sqrt(sum(Mod(c0)^2))

	eta <- matrix(0,nrow=M+1,ncol=1)

	for (k in 1:subintervals)
		eta <- eta + as.complex(frequencies[k]/(t(Conj(c0))%*%statisticsarray[,,k]%*%c0))*(statisticsarray[,,k]%*%c0)

	eta <- eta-c0

	newtonmanifold<-(c0+eta)
	newtonmanifold<-newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
	newtonmanifold<-newtonmanifold*exp(-1i*Arg(newtonmanifold[1]))
	newtonmanifoldprevious<-newtonmanifold

	for (j in 1:iter)
	{
	     eta<-matrix(0,nrow=M+1,ncol=1)
	     for (k in 1:subintervals)
             {
		eta <- eta + as.complex(frequencies[k]/(t(Conj(newtonmanifold))%*%statisticsarray[,,k]%*%newtonmanifold))*(statisticsarray[,,k]%*%newtonmanifold)
	     }
	     eta<-eta-newtonmanifold
	     #sqrt(sum(Mod(eta)^2))
	     newtonmanifold<-newtonmanifold+eta
	     #sqrt(sum(Mod(newtonmanifold)^2))
	     newtonmanifold<-newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
	     newtonmanifold<-newtonmanifold*exp(-1i*Arg(newtonmanifold[1]))
	     if (j==iter)
	          normsequence<-(sqrt(sum(Mod(newtonmanifold-newtonmanifoldprevious)^2)))
	     newtonmanifoldprevious<-newtonmanifold

	}

	newtonmanifold<-newtonmanifold/sqrt(2*pi)

	loglik <- nntsloglikInterval0to2pi(data,cutpoints,newtonmanifold,M) 
	#AIC <- -2*loglik + 2*(2*M)
	#BIC <- -2*loglik + (2*M)*log(sum(data))
	gradnormerror <- normsequence

	cestimatesarray <- data.frame(cbind(0:M,newtonmanifold))
	cestimatesarray[,1]<-as.integer(Re(as.matrix(cestimatesarray[,1])))
	names(cestimatesarray)[1]<-"k"
	names(cestimatesarray)[2]<-"cestimates"

	#res<-list(cestimates=cestimatesarray, loglik=loglik, AIC=AIC, BIC=BIC, gradnormerror=gradnormerror)
	res<-list(cestimates=cestimatesarray, loglik=loglik, gradnormerror=gradnormerror)
	return(res)
}

