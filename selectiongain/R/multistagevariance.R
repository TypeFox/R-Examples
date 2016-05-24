`multistagevariance` <-
function(Q,corr,alg=GenzBretz())
{
lim.y=-200
dim=length(Q)+1
k=c(lim.y,Q) 
alphaofx=pmvnorm(lower=k,corr=corr)

# in the case of one stage selection, dim = 2

for (i in 1:dim)
{
   if (is.infinite(k[i])==TRUE)
   {
     k[i]=-100
   }
}

# if dim < 2 , then stop

if (dim<2)
{
stop("dimension of k must bigger than 1, otherwise this is not one-stage selection")

}else if(dim==2)
{

# if dim=2, then we deal with one stage selection

gainresult=(1-corr[1,2]^2*(dnorm(k[2])/alphaofx))*(dnorm(k[2])/alphaofx- k[2]) 

}else
{

# otherwise, we use Talis(1965)'s algorithm 

A=array(0,c(dim,dim))

for (i in 1 : dim)
{  
   for (j in 1 : dim)
    {
          if(i!=j)
         {
              A[i,j]= (k[j]-corr[i,j]*k[i])/ (1-corr[i,j]^2)^0.5
          }
    }
}




part.corr=array(1,c(dim,dim,dim))

for (i in 1 : dim)
{  
   for (j in 1 : dim)
    {
         for (q in 1 : dim)
         {
         if(i!=j && q!=j && i!=q)

            { 
            part.corr[i,j,q]= (corr[i,j]-corr[i,q]*corr[j,q])/ ((1-corr[i,q]^2)^0.5 * (1-corr[j,q]^2)^0.5)
             }
          }

      }
}


BSR.Q=array(1,c(dim,dim,dim))

for (s in 1 : dim)
{  
   for (r in 1 : dim)
    {
         for (q in 1 : dim)
         {
         if(s!=r && q!=r && s!=q)

            { 
              BSR.Q[s,r,q]= (corr[r,s]-corr[r,q]*corr[s,q])/ (1-corr[r,q]^2)
             }
          }

      }
}


BSQ.R=array(1,c(dim,dim,dim))

for (s in 1 : dim)
{  
   for (r in 1 : dim)
    {
         for (q in 1 : dim)
         {
         if(s!=r && q!=r && s!=q)

            { 
              BSQ.R[s,q,r]= (corr[q,s]-corr[r,q]*corr[s,r])/ (1-corr[r,q]^2)
             }
          }

      }
}


RSR.Q=array(1,c(dim,dim,dim))

for (s in 1 : dim)
{  
   for (r in 1 : dim)
    {
         for (q in 1 : dim)
         {
         if(s!=r && q!=r && s!=q)

            { 
              
               RSR.Q [s,r,q]= (corr[s, r] - corr[s, q] * corr[q, r])/sqrt((1 - corr[s, q]^2) * (1 - corr[q, r]^2))
          
             }
          }

      }
}




A2=array(0,c(dim,dim,dim))

for (q in 1 : dim)
{  
   for (s in 1 : dim)
    {
      for (r in 1:dim)
        {
          if(q!=s && q!=r&& r!=s )
         {
              A2[q,s,r]= (k[s]-BSQ.R[s,q,r]*k[q]-BSR.Q[s,r,q]*k[r])/ ((1-corr[s,q]^2)*(1-RSR.Q[s,r,q]^2))^0.5
          }
        }  
    }
}




part.corr.second=array(1,c(dim,dim,dim,dim))

for (i in 1 : dim)
{  
   for (j in 1 : dim)
    {
         for (q in 1 : dim)
         {
            for (r in 1:dim)
             {
                    if(i!=j && q!=j && i!=q && r!=i&& r!=j&& r!=q)

                  { 
                  part.corr.second[i,j,q,r]= (part.corr[i,j,q]-part.corr[i,r,q]*part.corr[j,r,q])/ ((1-part.corr[i,r,q]^2)^0.5 * (1-part.corr[j,r,q]^2)^0.5)
                                }
             
             }
          }

      }
}




j3q<-function (q,A,part.corr,dim)

{      
    
       lower=A[q,-q]     
   
    corr= part.corr[-q,-q,q]
 
  
     
    output=pmvnorm(lower = lower, upper = rep(Inf,c(dim-1)), mean = rep(0, length(lower)), 
    corr = corr, sigma = NULL, algorithm =  alg) 
    
    
    output
 }


j3q.second<-function (q,r,A2,part.corr,dim)

{      
    
          
  
    corr= part.corr.second[c(-q,-r),c(-q,-r),q,r]

    lower=A2[q,c(-q,-r),r]
    
    if (dim>3)
    {
    output=pmvnorm(lower = lower, upper = rep(Inf,c(dim-2)), mean = rep(0, length(lower)), 
    corr = corr, sigma = NULL, algorithm =  alg) 
    }else
    {
    output=pnorm(q=lower, mean = 0, sd = 1, lower.tail =FALSE, log.p = FALSE)
    }
    
    output
 }





 
calculatx2<-function(i,j,A,A2,part.corr,part.corr.second,dim,corr,k,alpha3)
{ 
   sum1=0
   sum2=0
  
   for (q in 1 : dim)
 {
   sum1= sum1+ corr[i,q]*corr[j,q]*dnorm(k[q])*k[q]*j3q(q,A,part.corr,dim)/alpha3
 }
   
     for (q in 1 : dim)
    {    
       for (r in 1:dim)
         {
            if (q!=j && i!=q && r!=i&& r!=j&& r!=q )
            {
            

            
            sum2= sum2+ corr[q,i]*dmvnorm(x=c(k[q],k[r]),sigma=corr[c(q,r),c(q,r)])/alpha3 * j3q.second(q,r,A2,part.corr,dim) *(corr[r,j]-corr[q,r]*corr[q,j]/corr[q,q])
            
            
            }
         }
   }
  
  
   output =corr[i,j] + sum1 + sum2


  c(output,sum1,sum2,alpha3)
} 
  

gainresult<-calculatx2(i=1,j=1,A=A,A2=A2,part.corr=part.corr,part.corr.second=part.corr.second,dim=dim,corr=corr,k=k,alpha3=alphaofx)

gain=multistagegain(Q= Q, corr=corr,alg=alg)

# V = E(x^2)-E(x)^2

variance=gainresult[1]-gain^2

}



variance[1]



}

