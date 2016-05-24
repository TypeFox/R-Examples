"mainSimul2" <-
function(nbBed,nbPatient,disXi,disYi,toplot=FALSE)
{
   	Tab<-array(0,1);
   	Temp<-array(0,1);
   	indiceT<-1;
 
   	for (i in 1:nbBed)
   	{
           Temp<-simul2(nbPatient,disXi,disYi,toplot);
           if(length(Temp)>=1 && Temp[1] !=0)
           {
               Tab[indiceT:(indiceT+length(Temp)-1)]<-Temp;
               indiceT<-indiceT+length(Temp);
           }
   	}

   

   #calculus of R and C_R with the equation (*) only valid if disXi is an exponentiel distribution and disYi too
   
   	lambda <- disXi$param[[1]];
   	mu <- disYi$param[[1]];
  	R <- mu;

	CR<-lambda/R*(1-lambda/(R-lambda));
	  
	if(toplot)
	{
		print("value of lambda is")
		print(lambda)
		print("value of mu is")
		print(mu)
		print("value of R is")
		print(R)
		print("value of C_R is")
		print(CR)
	}
	
	return (list(CR=CR,R=R,T=Tab));


}

