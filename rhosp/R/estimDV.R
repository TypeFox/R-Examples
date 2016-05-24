"estimDV" <-
function(fileName,toplot=TRUE,header=TRUE,ks=FALSE)
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
    
    
    #calculus of R and C_R with the equation (*) only valid if disXi is an exponentiel distribution and disP depends on Xi
   
   	lambdaTilde <- as.numeric(DV(T)[1]);
   	
	muTilde <- as.numeric(DV(T)[2]);  

	R<- muTilde;
	CR<-lambdaTilde/muTilde*(1-lambdaTilde/(muTilde -lambdaTilde));

    if(toplot)
    {
	    print("DV")
	   	print("value of lambdaTilde is");
	    print(lambdaTilde)
	    print("value of  muTilde is")
	    print(muTilde)
		print("value of R is")
		print(R)
		print("value of C_R is")
		print(CR)
	}


    if(ks){    return (list(CR=CR,R=R,T=T,lambdaHat=lambdaTilde, muHat=muTilde, lambdaEmp=lambdaTilde,muEmp=muTilde,KolmogorvSmirnovTest=kstest));
    }
    else{    return (list(CR=CR,R=R,T=T,lambdaHat=lambdaTilde, muHat=muTilde, lambdaEmp=lambdaTilde,muEmp=muTilde));
    }
}

