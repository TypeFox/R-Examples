filter.F<-function(data,F.low=1,F.up=3)
{
   if(!is.data.frame(data) || !is.numeric(c(F.low,F.up)) || (F.low<1.0 || F.up>3.0) || F.low>F.up) 
      stop("invalid input parameter(s) specification: check data/F.low/F.up")
  
   data[data$F>=F.low & data$F<=F.up,]
}
