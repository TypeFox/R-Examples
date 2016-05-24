"estimParam" <-
function(fileName,toplot=TRUE,header=TRUE,ks=FALSE,DV=FALSE)
{
    echantillon<-read.table(fileName,header)
#     A data frame is a list of variables of the same length with unique
#     row names, given class '"data.frame"'.
#	
#	we suppose header is : no	sejour	evt
    
    nbPatient<-length(echantillon["no"][[1]]);
    sejour<-echantillon["sejour"][[1]]; #X_i
    stayCumsum<-cumsum(sejour); #"T_i"
    evt<-echantillon["evt"][[1]]; #Z_i
    
    if(ks) kstest<-ks.test(sejour,"pexp");
    
    
	T<-array(0,1);
    
    #calcul de T
    iFirstEmpT<-1;
    iOldT<-1;
    for (i in 1:nbPatient)
    {
 	    if(evt[i]==1)
 	    {
	 	    if(iFirstEmpT==1)
	 	    {
		 	    T[iFirstEmpT]<-stayCumsum[i];
		 	    iFirstEmpT<-iFirstEmpT+1;
		 	    iOldT<-i;
	 	    }
	 	    else
	 	    {
		 	    T[iFirstEmpT]<-stayCumsum[i]-stayCumsum[iOldT];
		 	    iFirstEmpT<-iFirstEmpT+1; 
		 	    iOldT<-i;
	 	    }
	    }	
    }

    if(DV)
    {
	    lambdaTilde <- as.numeric(DV(T)[1]);
   	
		muTilde <- as.numeric(DV(T)[2]); 
	}
	else
	{
		lambdaTilde <-NULL;
		muTilde <-NULL;
	} 
    
    aux<-function (mu,x,z)
    {
        if(z==0)
        {
            res<- - mu*x*(1-z) 
        }
        else
        {
            res <-  z*log(1-exp(-mu*x))
        }
	return(res)
    }

   
    lnL <- function(mu,X,Z)
    {
	    if(length(X)!=length(Z))
	    	stop("size error")
	    	
	    res <-0;
	    for(i in 1:length(X))
	    	res[i] <- aux(mu,X[i],Z[i])
	    	
	    return(sum(res))
    }
    
    muHat <- optimize(lnL,X=sejour,Z=evt,lower=0,upper=10000,maximum=TRUE)$maximum;

    lambdaHat <- (length(sejour)-1)/(sum(sejour)); #the ESBVM estimator of lambda
    
   	R<- muHat;
	CR<-lambdaHat/muHat*(1-lambdaHat/(muHat -lambdaHat));
 	
	if(toplot)
        {
	   	print("parametric estimation")
	    print("value of lambda is")
		print(lambdaHat)
		print("value of R is")
		print(R)
		print("value of C_R is")
		print(CR)
	}



    if(ks){    return (list(CR=CR,R=R,T=T,lambdaHat=lambdaHat, muHat=muHat, lambdaEmp=lambdaTilde,muEmp=muTilde,KolmogorvSmirnovTest=kstest));
    }
    else{    return (list(CR=CR,R=R,T=T,lambdaHat=lambdaHat, muHat=muHat, lambdaEmp=lambdaTilde,muEmp=muTilde));
    }
}

