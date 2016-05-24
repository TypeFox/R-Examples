`pomdev.extra` <-
function(object1, object2, eps=10^-30, nrange=1000, 
		fullmsd=FALSE,  plotting=FALSE,...) 
{
	if(!is.numeric(object1)|!is.numeric(object2))
        stop("objects must be numeric ")
 
  if(any(!is.finite(object1))|any(!is.finite(object2)))
        stop("objects must have finite values ")
 	
	r1<-range(object1,na.rm=TRUE)
	r<-extendrange(r1,f=0.01)
	bw1<-bw.nrd0(object1)
	
	r2<-range(object2)
	
  d1<-density(object1,bw=bw1,from=r[1],to=r[2],n=nrange)
	o1<-d1$y
	w <- o1 < eps
	if (any(w)) o1[w] <- eps

	if(fullmsd){
		z <- matrix(NA, nrow=5, ncol=1)
		rownames(z) <- c("POMDEV","overlap","KLdiv","MSD","CrossMSD")
	}else{
		z <- matrix(NA, nrow=4, ncol=1)
		rownames(z) <- c("POMDEV","overlap","KLdiv","MSD")
	}
	d2<-density(object2,bw=bw1,from=r[1],to=r[2],n=nrange)
	o2<-d2$y
	w <- o2 < eps
	if (any(w)) o2[w] <- eps
	pdf2 <- approxfun( d2$x, o2, yleft=eps, yright=eps)
	
	z[1,] <-  - 2 * sum ( log(pdf2(object1)) )       #"raw" POMDEV
  z[2,] <- ( r1[1]<=r2[2] & r1[2]>=r2[1] )		 #overlap
	z[3,] <- sum(o1 * (-log(o2) + log(o1)))			 #KLdiv
	z[4,] <- (mean(object1)-mean(object2))^2		 #MSD   
	
	if(fullmsd) 
	{
		if(length(object2)*length(object1)<10e6){
				z[5,] <- mean(sapply(object2, object1 , FUN="-")^2)	#full cross MSD
			}else{
				tot<-0
				for(ll in 1:length(object2))
						tot<-tot+sum(sapply( object1 ,object2[ll], FUN="-")^2)
				z[5,] <- tot / length(object2)*length(object1)
			}
	}
	if(plotting)
	{
		old.par <- par(no.readonly = TRUE)
		
		par(mfrow=c(2,2))
		
		plot(seq(from=r[1],to=r[2],length=nrange),o1,type="l",main="Field pattern density",xlab="Focus variable",ylab="Density")
		plot(seq(from=r[1],to=r[2],length=nrange),o2,type="l",main="Model pattern density on field range",xlab="Focus variable",ylab="Density")
		plot(density(as.vector(object2)),main="Model pattern density on free range",xlab="Focus variable")
		plot(0:1,0:1,type="n",axes=FALSE,frame.plot=FALSE,xlab="",ylab="",...)
		text(0.5,0.5,paste("Information about the density function: \n bw=",round(bw1,4),"\n lower bound=",round(r[1],3),"\n upper bound=",round(r[2],3),"\n n points=",nrange,"\n\n Pomdev score=",round(z[1,],3),sep=" "))
				
		on.exit(par(old.par))
	}
	return(list(field_data=object1,sim_data=object2,kernel_estimator=list(bw=bw1,eps=eps,range=r, nrange=nrange),result=z))
}

