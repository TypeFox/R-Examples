"estimNonParam" <-
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
    
    #calculus of R and C_R with the equation (*) only valid if disXi is an exponentiel distribution and disP depends on Xi
   
   	 
	
    lambdaHat <- (length(sejour)-1)/sum(sejour); #the ESBVM of lambda
	
    eqStar<-function(r)
    {
        res <- 0;
        nbZi0 <- 0;
        for(i in 1:length(sejour))
        {
            if(evt[i]==0)
	    	{
                    nbZi0 <- nbZi0 +1;
                    res <- res + exp(r*sejour[i]);
                }
        }
        return( abs(res - (length(sejour))^2/(nbZi0)) );
    }
    
    R <- optimize(eqStar,muTilde,lower=0,upper=lambdaHat)$minimum;

    CR<-lambdaHat/R*(1-lambdaHat/(R-lambdaHat));
    
    if(toplot)
    {
        print("non parametric estimation")
        print("value of lambda is");
        print(lambdaHat)
		print("value of R is")
		print(R)
		print("value of C_R is")
		print(CR)
    }


    if(ks){    return (list(CR=CR,R=R,T=T,lambdaHat=lambdaHat, muHat=NULL, lambdaEmp=lambdaTilde,muEmp=muTilde,KolmogorvSmirnovTest=kstest));
    }
    else{    return (list(CR=CR,R=R,T=T,lambdaHat=lambdaHat, muHat=NULL, lambdaEmp=lambdaTilde,muEmp=muTilde));
    }
}

