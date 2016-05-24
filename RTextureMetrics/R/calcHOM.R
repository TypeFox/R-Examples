calcHOM <-
function(rawmat)
{
   #create matrix for contrast weights
   size<-dim(rawmat)[1]
   mathomweights<-matrix(0, nrow=size, ncol=size)  
  
   for(i in 1:size)  #by row
   {
      for(a in 1:size)  #by col
      {
          mathomweights[i,a]<-1/(1+(a-i)^2)
      }#eo for a
 
   }#eo for i

   #return(mathomweights)  #getestet am 19.5.2009

   hom<-rawmat*mathomweights   
   return(sum(hom))
}
