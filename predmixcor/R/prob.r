##############################################################################
#calculating probability of low correlation for given r and integrated over phi,
#the integration is done using simpson's rule 
#this function will be applied inside the corrected sampling with slightly 
#change such that it searches the past values and saves new calculated values
prob.cor.sum = function(indexalpha,phi.range,N,r,n = 5,flag = FALSE)
{
	 
        f = function(x)
        {
            phi = x
	    prob.cor(indexalpha = indexalpha, N = N,
	             a = r*phi,b = r*(1-phi),flag = flag)		    
        }
        
	simpson(f,phi.range[1],phi.range[2],n)/
	(phi.range[2] - phi.range[1])
}

##############################################################################
#calculating probability of low correlation for given r and phi, which will be
#read into a and b when the function is called
#N----A vector of 4 components, 
#N[1]=NO of X=0 in 1st group
#N[2]=NO of X=1 in 1st group
#N[3]=NO of X=0 in 2nd group
#N[4]=NO of X=1 in 2nd group
#a,b---the parameters in the prior for theta
#indexalpha is a matrix calculated by function ``bound'' which gives out the
#bound of n11(the second and third columns) for each n01(the first column,
#details go to the function bound

prob.cor <- function(indexalpha, N, a, b, flag = FALSE)
{
        if(a == 0 || b == 0)
           return( 1 )
        else
        {

        if(length(N)!=4)
        stop("Please enter 4 values for N")
        
        n0.c1 = N[1]
        n1.c1 = N[2]
        n0.c2 = N[3]
        n1.c2 = N[4]
        n = sum(N)
        lp = -Inf
        
        
        log.prob.c1 = matrix(0,n0.c1+1,n1.c1+1)
        log.prob.c2 = matrix(0,n0.c2+1,n1.c2+1)
        
        #calculationg the tables of probabilities for C1
        for( k.c1 in 0:(n0.c1+n1.c1) )
        for( n11.c1 in max(0,k.c1-n0.c1):min(n1.c1,k.c1) )
        {
                n01.c1 = k.c1-n11.c1
                if(n11.c1 == max(0,k.c1-n0.c1) )
                   log.prob.c1[n01.c1+1,n11.c1+1]=
                        LOG.prob(n01.c1,n11.c1,n0.c1,n1.c1,a,b)
                else
                   log.prob.c1[n01.c1+1,n11.c1+1] = 
		         log.prob.c1[n01.c1+2,n11.c1] +
                         log(n01.c1+1) + log(n1.c1-n11.c1+1) -
			 log(n0.c1-n01.c1)-log( n11.c1 )
        }

         #calculationg the tables of probabilities for C2
         #first to check whether it is similar to log.prob.c1
         if(n0.c2 == n0.c1 && n1.c2 == n1.c1)
                log.prob.c2 = log.prob.c1
         else 
                if(n0.c2 == n1.c1 && n1.c2 == n0.c1)
                        log.prob.c2 = t ( log.prob.c1 )
                else
                #we have to calculate log.prob.c2
                for( k.c2 in 0:(n0.c2+n1.c2))
                for( n11.c2 in max(0,k.c2-n0.c2):min(n1.c2,k.c2) )
                {
                n01.c2 = k.c2-n11.c2
                if(n11.c2 == max(0,k.c2-n0.c2) )
                       log.prob.c2[n01.c2+1,n11.c2+1] =
                       LOG.prob(n01.c2,n11.c2,n0.c2,n1.c2,a,b)
                else
                       log.prob.c2[n01.c2+1,n11.c2+1] =
                       log.prob.c2[n01.c2+2,n11.c2] +
                       log ( n01.c2+1 )+log ( n1.c2-n11.c2+1 ) +
                        -log ( n0.c2-n01.c2 ) - log( n11.c2 )
                }

    #adding the probabilities of all possible pairs which have correlation
    #less than alpha

   if(!flag)
   {
        for(i in 2:nrow(indexalpha) )
        {   #0: max( (indexalpha[i,2]-1) ,0 ) ,
            #min( indexalpha[i,3] + 1 , n1.c1 + n1.c2 ) : ( n1.c1 + n1.c2 )
            n01 = indexalpha[i,1]
            
            for( n11 in  seq(0, indexalpha[i,2]-1) )
                lp = LOG.add.exp( lp ,
                                  prob.add ( n01,n11,n0.c1,n1.c1,n0.c2,
				             n1.c2,log.prob.c1,log.prob.c2
					   )
                                )
        }
	#print(exp(lp))
        return( 1 - 2*exp( lp ) )
   }
   else
   {
        for(i in 1:nrow(indexalpha) )
        {
            n01 = indexalpha[i,1]
            for(n11 in indexalpha[i,2]:indexalpha[i,3])
                lp = LOG.add.exp(lp ,
                         prob.add( n01,n11,n0.c1,n1.c1,n0.c2,
			    n1.c2,log.prob.c1,log.prob.c2
				 )
                                 )
        }

        return( exp( lp ) )     
  }
 }
}

#####################################################################    
#this function adds probabilities( stored in log form) in prob.c1 and prob.c2 
#which make up(n01,n11), then also stored in log form
prob.add = function(n01,n11,n0.c1,n1.c1,n0.c2,n1.c2,log.prob.c1,log.prob.c2)
{
	LOG.sum.exp (
		log.prob.c1 [ max(0,n01-n0.c2) : min(n01,n0.c1) + 1,
   		 max ( 0,n11-n1.c2 ) : min ( n11,n1.c1 )  + 1 ] +
    		log.prob.c2 [ n01 - max (0,n01-n0.c2 ) : min( n01,n0.c1) + 1,
		 n11- max ( 0,n11 - n1.c2 ) : min ( n11,n1.c1 )  + 1 ]         
                    )
}

######################################################################
#this function calculates the log-probability of n01 and n11 
#by Polya Urn Scheme in each group, which is used to calculate 
#``log.prob.c1'' and ``log.prob.c2'' in function ``prob.cor''
LOG.prob = function(n01,n11,n0,n1,a,b)
{
    if( n01>n0 || n11>n1 )
        stop("Wrong Arguments in 'Prob'")
	
    n = n0 + n1
    log ( choose ( n0,n01 ) ) + log ( choose ( n1,n11 ) ) +
    lgamma ( a+n01+n11 ) - lgamma ( a ) +
    lgamma ( b+n-n01-n11 )  - lgamma ( b ) +
    lgamma ( a+b ) - lgamma ( a+b+n )   

}

