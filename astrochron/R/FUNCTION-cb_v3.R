### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### cb: take input vectors, column bind, and return as data.frame 
###                          (SRM: May 2, 2012; May 24, 2013; June 7, 2013)
### 
###########################################################################

cb<- function (a,b)
{
  a <- data.frame(a)
  numvar <- ncol(a)
  if(numvar > 1)
    {
      out <- a[b[1]]    
      for(i in 2:length(b))
       {
         out <- cbind(out,a[b[i]])
       }
     }
  if(numvar == 1) out <- cbind(a,b)
  return( data.frame (out) )

### END function cb
}
