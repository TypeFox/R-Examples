prepRetOtp <- function(output,kmin,kmax,nsol,Waldval,algorithm=c("anneal","genetic","improve","eleaps"),Numprb=NULL,Optimal=NULL,tl=NULL)
{
	 algorithm <- match.arg(algorithm) 
	 if (algorithm!="eleaps")  {
        	valores <- matrix(nrow=nsol,ncol=length(kmin:kmax),output[[1]])
        	variaveis <- array(output[[2]],c(nsol,kmax,length(kmin:kmax)))
        	bestval <- output[[3]]
        	bestvar <- t(matrix(nrow=kmax,ncol=length(kmin:kmax),output[[4]]))
        	output <- list(subsets=variaveis,values=valores,bestvalues=bestval,bestsets=bestvar,call=output[[5]])
	 }
         dimnames(output$subsets)<-list(paste("Solution",1:nsol,sep=" "),paste("Var",1:kmax,sep="."),paste("Card",kmin:kmax,sep="."))
         dimnames(output$values)<-list(paste("Solution",1:nsol,sep=" "),paste("card.",kmin:kmax,sep=""))
         names(output$bestvalues)<-paste("Card",kmin:kmax,sep=".")
         dimnames(output$bestsets)<-list(paste("Card",kmin:kmax,sep="."),paste("Var",1:kmax,sep="."))

         ndim <- kmax-kmin+1
	 dimind <- 1:ndim
	 nosol <- dimind[sapply(dimind,function(dim) all(output$subsets[,,dim]==0))]

	 if (length(nosol)>0)  {

		wrnmsg1 <- paste(algorithm,"was not able to find any reliable non-collinear subset")

		if (algorithm=="eleaps") 

			wrnmsgEl <- paste("within the specified time limit of ",tl," seconds.\nTo search for reliable solutions either increase the value of the function argument 'timelimit',\n or try one of the available meta-heuristics (anneal, genetic or improve).\n")

		if (length(nosol)==ndim)  {
			if (algorithm=="eleaps") warning(paste(wrnmsg1,wrnmsgEl)) 
			else warning(paste(wrnmsg1,".\n",sep="")) 
			return(NULL)
		}
		else {
			if (length(nosol)==1)  wrnmsg2 <- paste("of dimensionality",kmin+nosol-1)
			else  wrnmsg2 <- paste("of dimensionalities",paste(kmin+nosol-1,collapse=" "))
			if (algorithm=="eleaps") warning(paste(wrnmsg1,wrnmsg2,"\n",wrnmsgEl)) 
			else warning(paste(wrnmsg1," ",wrnmsg2,".\n",sep="")) 

			output$subsets <- output$subsets[,,-nosol,drop=FALSE]  
			output$values <- output$values[,-nosol,drop=FALSE]  
			output$bestsets <- output$bestsets[-nosol,,drop=FALSE]  
			output$bestvalues <- output$bestvalues[-nosol]  
		}
	 } 
	for (dim in 1:length(output$bestvalues))  {
		trialind <- 1:nsol
		nullsols <- trialind[sapply(trialind,function(ind) all(output$subsets[ind,,dim]==0))] 
		if (length(nullsols)>0)  output$values[nullsols,dim] <- 0.  
 	 }
         if (!is.null(Waldval))  {
		criterion <- "WALD"
		criterio <- 8
		validvalues <- output$values[output$values > 0.]
		output$values[output$values > 0.] <- rep(Waldval,length(validvalues)) - validvalues * Waldval
		output$bestvalues <- rep(Waldval,ndim) - output$bestvalues * Waldval  
	 }

	 if (algorithm=="eleaps")  {
		if (Optimal == FALSE) 
	    		warning("eleaps was not able to complete an exact search within the specified time limit of ",tl," seconds, \n and the optimality of the returned solutions cannot be guaranteed.\n To search for better solutions either increase the value of the function argument 'timelimit',\n or try one of the available meta-heuristics (anneal, genetic or improve).\n") 
	}
	else  if (Numprb==TRUE)  
		warning("Because of numerical problems caused by strong multicolinearity\n some subsets were excluded from the analysis.\nYou can try to increase the number of subsets to be compared by reducing the value\nof the function argument 'tolval', but the numerical accuracy of results may be compromised.\n")

	 output  #  return(output)
}
