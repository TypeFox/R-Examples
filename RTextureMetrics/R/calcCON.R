calcCON <-
function(rawmat)
{
   #create matrix for contrast weights
   size<-dim(rawmat)[1]
   matconweights<-matrix(0, nrow=size, ncol=size)  
  
   for(i in 1:size)  #by row
   {
      for(a in 1:size)  #by col
      {
          matconweights[i,a]<-(a-i)^2
      }#eo for a
 
   }#eo for i

   con<-rawmat*matconweights   
   return(sum(con))
   
}
