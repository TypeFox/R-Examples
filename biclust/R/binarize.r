#Binarize by the input threshold or using median if no treshold given
binarize=function(x,threshold=NA)
  {
  matd=x
  if(is.na(threshold))
    {
    threshold=min(x)+(max(x)-min(x))/2
    print(paste("Threshold: ",threshold))
    }

  matd[matd<=threshold]=0
  matd[matd>threshold]=1
  matd
  }
  
#Binarize to get the input percentage of 1s over 0s
#The algorithm stops when it gets a density in [percentage-error, percentage+error]
#Gap is the increment taken in iterative search
binarizeByPercentage=function(x,percentage,error=0.2,gap=0.1)
  {
  if(error<=gap)
    {
    print("Error: error value must be greater than gap")
    }
  else
    {
    threshold=min(x)+(max(x)-min(x))/2
    dens=percentage
    repeat
      {
      matd=binarize(x,threshold)
      if(dens==densityOnes(matd))  {break}
      dens=densityOnes(matd)
      if(dens>percentage)
        {
        threshold=threshold+gap
        }
      else
        {
        threshold=threshold-gap
        }

      if(dens>=(percentage-error) && dens<=(percentage+error))
        break
      }
    print(paste("Threshold applied is ",threshold))
    matd
    }
  }

#---------------------- DENSITY -------------------------
#Percentage of number of 1s over number of 0s
densityOnes=function(x)
  {
  num1=length(x[x==1])
  num0=length(x[x==0])
  if((num1+num0)!=prod(dim(x))) print("Error: x must be a binary matrix")
  else
    {
    den=(num1/num0)*100  #Density of 1s over 0s
    den
    }
  }