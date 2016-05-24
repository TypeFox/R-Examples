##pattern cluster, return pattern as matrix
## equal to the function:
## unique(X.discret)
patternmatrix<-function(unipattern,n.entity)
{
  len <- n.entity+1
  px<-matrix(0,nrow=length(unipattern),ncol=n.entity)

  for(i in 1:length(unipattern))
  {
     n<-nchar(unipattern[i])
      for(j in 1:n)
      {
       px[i,len-j]<-as.integer(substr(unipattern[i],n-j+1,n-j+1))
      }
   }

   rownames(px)<-unipattern
   return(px)
 }

