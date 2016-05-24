"DV" <-
function(T)
{
   ksi1 <- mom1(T);
   ksi2 <- mom2(T);

   	delta <- 2*ksi2-3*(ksi1^2);
    if (delta > 0)
    {
    	
   		M <- min( ((-1)*ksi1 + sqrt(abs(delta)))/(ksi2-2*(ksi1^2)) , ((-1)*ksi1 - sqrt(abs(delta)))/(ksi2-2*(ksi1^2)) );
   		if(M<0)	M<- max( ((-1)*ksi1 + sqrt(abs(delta)))/(ksi2-2*(ksi1^2)) , ((-1)*ksi1 - sqrt(abs(delta)))/(ksi2-2*(ksi1^2)) );
   		if(M<0)	stop("error: root M negative")
    	muEmp <- M  
   		lambdaEmp <- muEmp/(muEmp*ksi1-1);   
    }
    else 
    {
        stop("error:delta negative");
    }
    
	
   return (c(lambdaEmp,muEmp));

}

