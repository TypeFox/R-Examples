calcDIS <-
function(rawmat)
{
   #create matrix for contrast weights
   size<-dim(rawmat)[1]
   matconweights<-matrix(0, nrow=size, ncol=size)  
  
   for(i in 1:size)  #by row
   {
      for(a in 1:size)  #by col
      {
          matconweights[i,a]<-abs(a-i)
      }#eo for a
 
   }#eo for i

   dis<-rawmat*matconweights   
   return(sum(dis))
}
