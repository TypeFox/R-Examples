"calcErrorDV" <-
function(file,nb=10,disXi,disP,plot=TRUE)
{
    tabREstim <- 0;
    tabCREstim <- 0;

    tabRTh <- 0;
    tabCRTh <- 0;
    

    for (i in 1:nb)
    {
	    makeSample(file,1000,disXi,disP)
        
        #calculus of R and C_R with the equation (*) only valid if disXi is an exponentiel distribution and disP depends on Xi
	   
	   	if(disP$nbparam==1) q<-function(x){1-disP$disfun(x,disP$param[[1]])};
	    if(disP$nbparam==2) q<-function(x){1-disP$disfun(x,disP$param[[1]],disP$param[[2]])};
	    if(disP$nbparam==3) q<-function(x){1-disP$disfun(x,disP$param[[1]],disP$param[[2]],disP$param[[3]])};
	
	   	lambda <- disXi$param[[1]]
	  	Lq <- function(s)
	  	{
	       	g <- function(y,t)	q(y)*exp(-y*t)    
	          
	       	return(integrate(g,0,+Inf,t=s)$value);	
	  	}
	    
	  	Eq <- function(r){	abs(lambda*Lq(lambda-r)-1);	}
 		tabRTh[i] <- optimize(Eq,0,lower=0,upper=lambda)$minimum;
		tabCRTh[i] <- lambda/tabRTh[i]*(1-lambda/(tabRTh[i]-lambda));
        
        ei<-estimDV(file,toplot=FALSE,header=TRUE)
        tabREstim[i] <- ei$R;
        tabCREstim[i] <- ei$CR;
    }
    if(plot)
    {
        print("bias is");
        print(mean(tabREstim-tabRTh));
        print("variance is");
        print(sd(tabREstim)^2);
        print("mean value of R is");
        print(mean(tabREstim));
        print("mean value of C_R is");
        print(mean(tabCREstim));
    }
    
   
return(list(bias=mean(tabREstim-tabRTh),var=sd(tabREstim)^2,R=mean(tabREstim),CR=mean(tabCREstim)))

}

