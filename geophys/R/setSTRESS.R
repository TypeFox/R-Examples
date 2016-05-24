setSTRESS<-function(Stensor)
{

  
  if(is.matrix(Stensor))
    {
      ES = eigen(Stensor)
    }

  
  if(is.vector(Stensor))
    {
      if(length(Stensor)==3)
        {
          Stensor = diag(Stensor)  
        }

      
      if(length(Stensor)==6)
        {
          nstress = Stensor[1:6]
          MAT = matrix(ncol=3, nrow=3)

          MAT[ lower.tri(MAT, diag=TRUE)] = nstress
          
          MAT[ upper.tri(MAT) ] = MAT[lower.tri(MAT)]

          Stensor = MAT


        }
      

      if(any(is.na(Stensor)))
        {
          print("BAD Stress Tensor!")
          return(NULL)
        }
      else
        {
          ES = eigen(Stensor)
        }
    }

  return(ES)

}
