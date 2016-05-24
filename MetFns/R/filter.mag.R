filter.mag<-function(data,mag.low=3,mag.up=8)
{
   if(!is.data.frame(data) || !is.numeric(c(mag.low,mag.up)) || (mag.low<3 || mag.up>8) ||  mag.low>mag.up) 
      stop("invalid input parameter(s) specification: check data/mag.low/mag.up")
  
   data[data$lmg>=mag.low & data$lmg<=mag.up,]
}
