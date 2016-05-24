# mutistagecor.R author: xuefei mi, 11-03-2013, for selectiongain package v2.0.2

# modified at 2015-11-09, funtion is modifyed to use multiple cpu cores


`multistagegain` <-
function(corr, Q, alg= GenzBretz(),parallel=FALSE)
{
  
  if (parallel)
  {    
  #library(parallel)
  
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores);
#  clusterEvalQ(cl,library(mvtnorm))
  }
  
# internal default parameters
  Vg=1
  partial=FALSE
  stages=partial
  lim.y=-200
  k=c(lim.y,Q)
  sum.dim=length(k)
  alphaofx=pmvnorm(lower=k,corr=corr,algorithm=alg)
  dim=sum.dim

# main function begins

# check if Q and corr have the same dimension

  if (length(Q)==dim(corr)[1]-1)
  {
  }else if (length(Q)!=dim(corr)[1]-1) 
  {
    stop("dimension of Q must be same as dim(corr)[1]+1")
  }

# check if k[i]= inf

  for (i in 1:dim)
  {
     if (is.infinite(k[i])==TRUE)
     {
       k[i]=-lim.y
     }
  }

# check dim

  if (dim<2)
  {
    stop("dimension of k must bigger than 1, otherwise this is not one-stage selection")
  }else if(dim==2)
  {

# in the case of one stage selection, dim = 2

    gainresult=corr[1,2]*dnorm(k[2])/alphaofx
  }else
  {

# calculate the gain according to Talis(1965)

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

      j3q<-function (q,A,part.corr,dim,alg)
      {            
        lower=A[q,-q]
        corr= part.corr[-q,-q,q]
        output=pmvnorm(lower = lower, upper = rep(Inf,c(dim-1)), mean = rep(0, length(lower)),   corr = corr, sigma = NULL, algorithm =  alg) 
        output
       }

      loopi<-function(i,k,A,part.corr,dim,alg,alpha3)
      {
      outputarrayi<-corr[1,i]*dnorm(k[i])*j3q(i,A,part.corr,dim,alg)/alpha3
      outputarrayi
      }
      
    
      
      
      calculatx1<-function(A,part.corr,dim,corr,k,alpha3,stages=FALSE)
      {  
        if (stages)
        {
           outputarray=rep(0,dim)
          i=1
          
          if (parallel)
          {
            outputarray <- parSapply(cl=cl, 1:dim, FUN=loopi,k=k,A=A,part.corr=part.corr,dim=dim,alg=alg,alpha3=alpha3)
            outputarray[2:dim]
          }else
          { 
          for (i in 1 : dim)
          {
            outputarray[i]=corr[1,i]*dnorm(k[i])*j3q(i,A,part.corr,dim,alg)/alpha3
           }
          outputarray[2:dim]
            }  
      
        }else
        {   
         outputarray=rep(0,dim)
          if (parallel)
          { 
            
            #outputarray <- parLapply(cl=cl, 1:dim, fun=loopi,k,A,part.corr,dim,alg,alpha3)
            
            outputarray <- parSapply(cl=cl, 1:dim, FUN=loopi,k=k,A=A,part.corr=part.corr,dim=dim,alg=alg,alpha3=alpha3)
            
            output<-sum(outputarray) 
            output
          }else
          {
          
            output=0
            i=1
          
         
          
           for (i in 1 : dim)
           {
             outputarray[i]=corr[1,i]*dnorm(k[i])*j3q(i,A,part.corr,dim,alg)/alpha3
           }
            sum(outputarray)
          }
         
        }
      }
      gainresult<-calculatx1(A=A,part.corr=part.corr,dim=dim,corr=corr,k=k,alpha3=alphaofx,stages=stages)
    }
  if (parallel)
  { 
  stopCluster(cl);
  }
  if (stages==TRUE)
  {
    gainresult*Vg^0.5
  }else
  {
    gainresult[[1]]*Vg^0.5
  }

}

