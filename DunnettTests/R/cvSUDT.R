cvSUDT <-
function(k,alpha=0.05,alternative="U",df=Inf,corr=0.5,corr.matrix=NA,mcs=1e+05){
  if(k > 16){
  	stop("The funtion is not applicable to the situations where the number of tests exceeds 16")
  }
  
  cvSet <- NULL
  z0 <- rnorm(mcs)
  if(df == Inf){u <- 1}else{u <- sqrt(rchisq(mcs,df=df)/df)}
  
  #for different tail test, the function for calculate Psi(di)
  if(alternative=="U"){
  	psi.fun <- function(cj,corr){pnorm((cj*u+z0*sqrt(corr))/sqrt(1-corr))}
    }
  if(alternative=="B"){
    psi.fun <- function(cj,corr){pnorm((cj*u+z0*sqrt(corr))/sqrt(1-corr))-pnorm((-cj*u+z0*sqrt(corr))/sqrt(1-corr))}
    }
  
  ##calculate c1
  if(alternative=="U"){
    c1 <- qt(1-alpha,df)
  }
  if(alternative=="B"){
    c1 <- qt(1-alpha/2,df)
  }
  cvSet <- c(cvSet,c1)
  J0 <- 1
  
  #for equal correlations
  if(is.matrix(corr.matrix)==FALSE){
    #the function for finding cj for equal corr
    cj.fun.balanced <- function(cj,corr,J.fun,list.J,list.Psi){
  	 	Psi.dj <- psi.fun(cj,corr)
  	 	mean(J.fun(Psi.dj,list.J,list.Psi))-(1-alpha)
  	 	}    
  	
  	##calculate c2
    Psi.d1 <- psi.fun(c1,corr)
  	J1 <- Psi.d1
  	list.J <- c(list(J0),list(J1))
  	list.Psi <- list(Psi.d1)	
  	c2 <- uniroot(cj.fun.balanced,corr=corr,J.fun=J2.fun,list.J=list.J,list.Psi=list.Psi,interval= c(0,10),tol=1e-05)$root 
  	cvSet <- c(cvSet,c2)
 	
  	##calculate c3-cvSet
  	if(k >= 3){
  	list.J.Fun <- list(J2.fun,J3.fun,J4.fun,J5.fun,J6.fun,J7.fun,J8.fun,J9.fun,J10.fun,J11.fun,J12.fun,J13.fun,J14.fun,J15.fun,J16.fun)
  	for(t in 3:k){
  	    Psi.dts <- psi.fun(cvSet[t-1],corr)
        Jts <- list.J.Fun[[t-2]](Psi.dts,list.J,list.Psi)  
        list.J <- c(list.J,list(Jts))
        list.Psi <- c(list.Psi,list(Psi.dts))
        ct <- uniroot(cj.fun.balanced,corr=corr,J.fun=list.J.Fun[[t-1]],list.J=list.J,list.Psi=list.Psi,interval= c(0,10),tol=1e-05)$root
  		cvSet <- c(cvSet,ct)
  		}    
  	}
  	}

  #when the correlations are not equal
  if(is.matrix(corr.matrix)){
  	corr.min <- NULL
  	corr.max <- NULL
  	for(i in 2:k){
  		#find all the subsets of test statistics with subset size be i
  		subset <- combn(1:k,i)
  		#the averaged correlation for each subset of test statistics
  		ave.rho <- apply(subset,2,function(x){
  			corr.subset <- corr.matrix[x,x]
  			sum(corr.subset[row(corr.subset)<col(corr.subset)])/(i*(i-1)/2)
  		})
        #find the minimum and maximum of ave.rho over all the subsets
  		corr.min[i-1] <- min(ave.rho)
  		corr.max[i-1] <- max(ave.rho)
        }
    cj.fun.unbalanced <- function(cj,corr.min,corr.max,J.fun,list.J.min,list.Psi.min,list.J.max,list.Psi.max){
    	Psi.dj.min <- psi.fun(cj,corr.min)
        prob.min <- mean(J.fun(Psi.dj.min,list.J.min,list.Psi.min))
        Psi.dj.max <- psi.fun(cj,corr.max)
        prob.max <- mean(J.fun(Psi.dj.max,list.J.max,list.Psi.max))
        min(prob.min, prob.max)-(1-alpha)
    }
    ##calculate c2
    Psi.d1.min <- psi.fun(c1,corr.min[1])
  	J1.min <- Psi.d1.min
  	list.J.min <- list(J0,J1.min)
  	list.Psi.min <- list(Psi.d1.min)
  		
  	Psi.d1.max <- psi.fun(c1,corr.max[1])
  	J1.max <- Psi.d1.max
  	list.J.max <- list(J0,J1.min)
  	list.Psi.max <- list(Psi.d1.max)
     
  	c2 <- uniroot(cj.fun.unbalanced,corr.min=corr.min[1],corr.max=corr.max[1],
  	J.fun=J2.fun,list.J.min=list.J.min,list.Psi.min=list.Psi.min,
  	list.J.max=list.J.max,list.Psi.max=list.Psi.max,interval= c(0,10),tol=1e-05)$root
  	cvSet <- c(cvSet,c2)
  	
  	##calculate c3-cvSet
  	if(k >= 3){
  	list.J.Fun <- list(J2.fun,J3.fun,J4.fun,J5.fun,J6.fun,J7.fun,J8.fun,J9.fun,J10.fun,J11.fun,J12.fun,J13.fun,J14.fun,J15.fun,J16.fun)
  	for(t in 3:k){
  	    #Notes: for the calculation of ct, t=3,...,k, all of Psi(di) and Ji should be recalculated for the new corr.min[t-1] and corr.max[t-1]
  	    Psi.d1.min <- psi.fun(c1,corr.min[t-1])
        J1.min <- Psi.d1.min
        list.Psi.min <- list(Psi.d1.min)
        list.J.min <- list(J0,J1.min) 
        Psi.d1.max <- psi.fun(c1,corr.max[t-1])
        J1.max <- Psi.d1.max
        list.Psi.max <- list(Psi.d1.max)
        list.J.max <- list(J0,J1.max)
        for(i in 2:(t-1)){
        	Ji.fun <- list.J.Fun[[i-1]]
        	Psi.di.min <- psi.fun(cvSet[i],corr.min[t-1])
        	Ji.min <- Ji.fun(Psi.di.min,list.J.min,list.Psi.min)
            list.Psi.min <- c(list.Psi.min,list(Psi.di.min))
        	list.J.min <- c(list.J.min,list(Ji.min))
        	
        	Psi.di.max <- psi.fun(cvSet[i],corr.max[t-1])
        	Ji.max <- Ji.fun(Psi.di.max,list.J.max,list.Psi.max)
            list.Psi.max <- c(list.Psi.max,list(Psi.di.max))
        	list.J.max <- c(list.J.max,list(Ji.max))
        }
  		ct <-uniroot(cj.fun.unbalanced,corr.min=corr.min[t-1],corr.max=corr.max[t-1],J.fun=list.J.Fun[[t-1]],list.J.min=list.J.min,list.Psi.min=list.Psi.min,list.J.max=list.J.max,list.Psi.max=list.Psi.max,interval=c(0,10),tol=1e-05)$root
  		cvSet <- c(cvSet,ct)
  	}
}
}
names(cvSet)=paste("c",seq(1,k,by=1),sep="")
return(round(cvSet,3))
}
