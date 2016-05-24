filter.sol<-function(data, sol.low=0, sol.up=360)
{
   if(!is.data.frame(data) || !is.numeric(c(sol.low,sol.up)) || any(c(sol.low,sol.up)<0) || any(c(sol.low,sol.up)>360)) 
      stop("invalid input parameter(s) specification: check data/sol.low/sol.up ")

   year<-data$year[100]
   year<-ifelse(year<20,year+2000,1900+year)

   if(sol.up>=sol.low)
     data[data$sollong>=sol.low & data$sollong<=sol.up,]
   else if(sol.low>sol.up&&sol.up<=solar.long(year,12,31,23.99)&&sol.low>=solar.long(year,1,1,0)) 
     rbind(data[data$sollong>=sol.low,],data[data$sollong<=sol.up,])
   else
     stop("invalid input parameter(s) specification: check sol.low/sol.up format")
}
