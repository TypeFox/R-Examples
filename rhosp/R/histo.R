"histo" <-
function(X,disXi=NULL,disP=NULL,plotDV=FALSE)
{      
       arg2Cst<-list(disfun=pcst<-function(x,p) p ,nbparam=1,param=list(1/20));
    
	if(	length(X$T) !=0)
	{
		#initialise variables
	    C <- X$CR;
	    R <- X$R;
	    x <- X$T;
	    lambdaEmp <- as.numeric(X$lambdaEmp);
	    muEmp <- as.numeric(X$muEmp);
        lambdaHat <- as.numeric(X$lambdaHat);
        muHat <- as.numeric(X$muHat);
		
		
	    y <- sort(x);
        n <- length(x);
        k <- round(1+ log(n)/log(2)); # number of rectangles
        size <- n/k ; # size: number of points per rectangle

        #bounds of histogram
        bounds <- array (0, k+1);
        bounds[1] <- y[1] - .025*(y[n]-y[1]);
        bounds[k+1] <- y[n] + .025*(y[n]-y[1]);
        
        for(i in 2:k)   bounds[i] <- (y[size*(i-1)] + y[size*(i-1)+1])/2;	#class of the same size

        # plot the histogram
        
        donnee<-hist(x, prob=TRUE, breaks = bounds,  xlim=c(bounds[1], bounds[k+1]) );    
       	
		# plot the empiric density
       	#lines(density(x,kernel="epanechnikov"),col="blue")

       	# density of the random variable associated to the equivalent and not the distribution tail
   		densityEquivalent <- function(y) R*C*exp(-R*y);	
        lines(seq(0,max(x),by=0.1),densityEquivalent(seq(0,max(x),by=0.1)),col="green",lwd=1);
        
        # theoric density, survival function and relative distance
     
        thDenDef<-FALSE;
        funSurDef<-FALSE;
       
        if(!identical(disXi,NULL) && !identical(disP,NULL) && identical(disXi$rangen,rexp) && identical(disP$disfun,pexp)&& !identical(lambdaHat,numeric(0)) && !identical(muHat,numeric(0)))
        {
	       	survivalFun <- function(x)
	       	{
			        lambda <- lambdaHat;
			        mu <- muHat;
			   		-1/(lambda-mu)*(lambda *exp(-mu*x)-mu*exp(-lambda*x));	       
		   	}
                funSurDef <- TRUE;
                lines(seq(0,max(x),by=0.1),survivalFun(seq(0,max(x),by=0.1)),col="blue");
        }
		#print(identical(disXi$rangen,rexp))
		#print(identical(disP$disfun,pexp))
		#print(disP$disfun)
		#print(pexp)
        if(!identical(disXi,NULL) && !identical(disP,NULL) && identical(disXi$rangen,rexp) && identical(disP$disfun,pexp))
        {
	       	densityTh <- function(x)
	       	{
			        lambda <- disXi$param[[1]];
			        mu <- disP$param[[1]];
			   		-lambda*mu/(lambda-mu)*(exp(-lambda*x)-exp(-mu*x));	       
		   	}
                thDenDef<-TRUE;
        }
        
        if(identical(disXi$rangen,rexp) && identical(disP$disfun,arg2Cst$disfun))
		{
	    	
    	    	# case where xi has an exponential distribution and p is a constant
	    	densityTh <- function(x)
	       	{
			        
	       	    lambda <- disXi$param[[1]];
	            p <- disP$param[[1]];
		    	lambda*p*exp(-lambda*p*x);	       
	   		}
                thDenDef<-TRUE;
        }
        if(plotDV && !identical(lambdaEmp,numeric(0)) && !identical(muEmp,numeric(0)))
        {
                densityDV <- function(x)
				{
        	        return(-lambdaEmp*muEmp/(lambdaEmp-muEmp)*(exp(-lambdaEmp*x)-exp(-muEmp*x)));	       
	   			}
        }


        if(plotDV && thDenDef)
        {
            lines(seq(0,max(x),by=0.1),densityTh(seq(0,max(x),by=0.1)),col="red");
            lines(seq(0,max(x),by=0.1),densityDV(seq(0,max(x),by=0.1)),col="pink");
            #leg.txt <- c("Empiric density", "Equivalent density","Theoritical density","De Vielder's density")
            #legend("topright", legend = leg.txt, pch = 21,col=c(4,3,2,6))
            leg.txt <- c("Equivalent density","Theoritical density","De Vielder's density")
	    	legend("topright", legend = leg.txt, pch = 21,col=c(3,2,6))
	    	
        }
        else
        {
            if(plotDV)
            {
                lines(seq(0,max(x),by=0.1),densityDV(seq(0,max(x),by=0.1)),col="pink");
                #leg.txt <- c("Empiric density", "Equivalent density","De Vielder's density")
	      		#legend("topright", legend = c("Empiric density", "Equivalent density","De Vielder's density", pch = 21,col=c(4,3,2)))
                leg.txt <- c( "Equivalent density","De Vielder's density")
	      		legend("topright", legend = c("Empiric density", "Equivalent density","De Vielder's density", pch = 21,col=c(3,2)))
            }
            else
            {	
	            if(thDenDef)
		        {
		            lines(seq(0,max(x),by=0.1),densityTh(seq(0,max(x),by=0.1)),col="red");
		            
		            #leg.txt <- c("Empiric density", "Equivalent density","Theoritical density")
		            #legend("topright", legend = leg.txt, pch = 21,col=c(4,3,2,6))
		            leg.txt <- c("Equivalent density","Theoritical density")
			    	legend("topright", legend = leg.txt, pch = 21,col=c(3,2))
			    	
			    	
		        }
	            else
	            {
	                #leg.txt <- c("Empiric density", "Equivalent density")
		      		#legend("topright", legend = leg.txt, pch = 21,col=c(4,3))
	                leg.txt <- c("Equivalent density")
		      		legend("topright", legend = leg.txt, pch = 21,col=c(3))
      			}
            }

        }      
	}
	else
	{
		stop("empty data T")
	}
    return(NULL)
}

