"mainSimul" <-
function(nbBed,nbPatient,disXi,disP,toplot=FALSE,calc=TRUE)
{
   	Tab<-array(0,1);
   	Temp<-array(0,1);
   	indiceT<-1;
   	for (i in 1:nbBed)
   	{
           Temp<-simul(nbPatient,disXi,disP,toplot);
           if(length(Temp)>=1 && Temp[1] !=0)
           {
                Tab[indiceT:(indiceT+length(Temp)-1)]<-Temp;
               	indiceT<-indiceT+length(Temp);
           }
   	}


   if(calc	)
   {
	   #calculus of R and C_R with the equation (*) only valid if disXi is an exponentiel distribution and disP depends on Xi
	   
	   	if(disP$nbparam==1) q<-function(x){1-disP$disfun(x,disP$param[[1]])};
	    if(disP$nbparam==2) q<-function(x){1-disP$disfun(x,disP$param[[1]],disP$param[[2]])};
	    if(disP$nbparam==3) q<-function(x){1-disP$disfun(x,disP$param[[1]],disP$param[[2]],disP$param[[3]])};

	   	lambda <- disXi$param[[1]]
	  	Lq <- function(s)
	  	{
	       	g <- function(y,t)	q(y)*exp(-y*t)
	
               	# Class "integrate" with components from function integrate
				# value                   the final estimate of the integral.
				# abs.error               estimate of the modulus of the absolute error.
				# subdivisions    		the number of subintervals produced in the subdivision process.
				# message                 "OK" or a character string giving the error message.
				# call                    the matched call.	    
	          
	       	return(integrate(g,0,+Inf,t=s)$value);	
	  	}
	    
	  	Eq <- function(r)
	  	{
	           abs(lambda*Lq(lambda-r)-1);	# minimise la valeur absolue de la difference
		}	      		
	   
 		R <- optimize(Eq,1,lower=0,upper=lambda)$minimum;

		CR<-lambda/R*(1-lambda/(R-lambda));
	}
    else
    {
	 	R<-NULL
	 	CR<-NULL   
    }
      
	if(toplot && calc)
	{
		print("value of lambda is")
		print(lambda)
		if(isTRUE(identical(disXi$rangen,rexp)) && isTRUE(identical(disP$disfun,pexp)))
             {
			mu <- disP$param[[1]];
			print("value of mu is")
			print(mu)				
		}	
		print("value of R is")
		print(R)
		print("value of C_R is")
		print(CR)
	}
	
    #lambdaTilde <- as.numeric(DV(Tab)[1]);
   	lambdaTilde <-1
    #muTilde <- as.numeric(DV(Tab)[2]); 
	muTilde <-1
	
	return (list(CR=CR,R=R,T=Tab,lambdaEmp=lambdaTilde,muEmp=muTilde,lambdaHat=NULL,muHat=NULL));
}

