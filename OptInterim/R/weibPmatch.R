weibPmatch<-function(x,p0,shape,scale){
  if(missing(shape)& missing(scale))
    stop('One of shape or scale must be specified')
  if(!missing(shape)& !missing(scale))
    stop('One of shape or scale must be specified')

  nlogp0<- -log(p0)
  
  if(missing(shape)){
      return( log(nlogp0)/log(x/scale) )
  }else{
      return( x/nlogp0^(1/shape) )
  }
}
