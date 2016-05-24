

`multistageoptimum.nlm`<-function(corr, Vg=1, ini.value=c(N.upper+N.lower)/4, Budget, CostProd, CostTest, Nf=1, iterlim=100, alg=GenzBretz(),N.upper, N.lower)


{
   CostInitial=CostProd
   Vgca=Vg
   CostTv=CostTest  
   percent =0.0001
   N.fs=Nf

# main function begins
  if (length(CostTest)!= dim(corr)[1]-1)
  {
    stop( "dimension of CostTest has to be dim(corr)[1]-1")
  }  
  if (length(CostProd)!= dim(corr)[1]-1)
  {
    stop( "dimension of CostProd has to be dim(corr)[1]-1")
  }  


  if (Budget> sum( N.upper*(CostTest+CostProd)))
  {
    stop("N.upper is too small, try to set sum( N.upper*(CostTest+CostProd))>Budget")
  }
  	

  if (length(N.upper)>5 || length(N.upper)< 2)
  {
    stop( "dimension of N.upper should be between 2 and 5")
  }
   
  if (length(N.upper)!= length(N.lower))
  {
    stop( "dimension of upper and lower limit has to be equal")
  }
 
   if (length(N.upper)!= dim(corr)[1]-1)
  {
    stop( "dimension of upper has to be dim(corr)[1]-1")
  }

  
  
# check, if initial value is greater than budget

   cost = sum(ini.value*CostInitial)+sum(CostTv*ini.value)
   if (cost>Budget)
   {
    stop( "initial cost is higher than the budget, try other values")
   }

# MultiStageOptim() calls multistagegain() to calculate the selection gain for multi-stage selection
# here xN is a variable, and all the other parameters are fixed
# this is just a shell function, which fits the requirement of the optimization function
      
   MultiStageOptim<-function(xN,Vg,corr,N.fs,alg)       
{       
 output=0
 
 length.xN=length(xN)
      
 xN.matrix=embed(c(xN,N.fs),2)

    if (all(xN.matrix[,1] <= xN.matrix[,2]))
    {   
     alpha = xN.matrix[,1]/ xN.matrix[,2]
  
     Quantile= multistagetp(alpha, corr=corr, alg=alg)     
                           
     output<- -multistagegain(Q=Quantile, corr=corr, alg=alg)*Vg^0.5
     
    # if (cost<=Budget*0.9)
    # {
    #   output=output*0.8
    # }
     
     
    }#else
     #   { 
     #    output=0
     #   }

 

 output*10000
}
     
     
     
     # nlm algorithm
   
   matrix.a=array(0,c(length(CostProd),length(CostProd)))
     
     for (i in 1:length(CostProd))
     {
       for (j in 1:length(CostProd))
       {
         if (i-j==1)
         matrix.a[i,j]=-1
     
       }
     
     }
     
    matrix.b= diag(length(CostProd)) + matrix.a
     
 Amatrix= t(cbind(matrix.b,-CostProd-CostTest))

 
 Bvec=c(rep(0,length(CostProd)),-Budget)    
     
                                       
result=constrOptim(theta=ini.value,grad=NULL, f=MultiStageOptim, ui=Amatrix, ci=Bvec, Vg=Vg, corr=corr,N.fs=N.fs,alg=alg)
  
result.table = array(0,length(N.upper)+1)

sample.size=result$par

result.table=c(sample.size, -result$value/10000)
 
names(result.table)=c(rownames(sample.size, do.NULL = FALSE, prefix = "N"),"value")

result.table
  
}

