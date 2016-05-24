"keyfct"<-function(distance, key.scale,key.shape, ftype)
# Simple function to switch between the key functions
{

#   if(ftype=="hn"){
      res<-keyfct.hn(distance,key.scale)
#   }else if(ftype=="hz"){
#      res<-keyfct.hz(distance,key.scale,key.shape)
#   }

   return(res)

}
