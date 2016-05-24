filter<-function(data,year=NULL,month.beg=NULL, month.end=month.beg,day.beg=NULL,day.end=day.beg,time.beg=0,time.end=2359,
shw=NULL, imocode=NULL,long.low=0,long.up=180,ew=c("E","W"), lat.low=0,lat.up=90,ns=c("N","S"),name=NULL,fname=NULL, 
site=NULL,country=NULL,mag.low=3, mag.up=8,F.low=1.0, F.up=3.0,sol.low=0,sol.up=360,Ralpha=NULL,Delta=NULL,
h.low=10,h.up=90,r=NULL, C=5)
{
   data.select<-data

   if(!is.null(shw))
      data.select<-filter.shw(data.select, shw)   
   
   if(!is.null(year) && !is.null(month.beg) && !is.null(day.beg))
      data.select<-filter.datetime(data.select,year,month.beg, month.end, day.beg,day.end,
                                   time.beg, time.end)
 
   if(!is.null(imocode)) 
      data.select<-filter.imocode(data.select,imocode)
 
   if(!is.null(name) && !is.null(fname)) 
      data.select<-filter.obsname(data.select,name,fname)
 
   if(any(c(long.low,lat.low)>0) || long.up!=180 || lat.up!=90)
      data.select<-filter.gc(data.select,long.low,long.up,ew,lat.low,lat.up,ns)
  
   if(!is.null(site)) 
      data.select<-filter.site(data.select,site)
 
   if(!is.null(country)) 
      data.select<-filter.country(data.select,country)
 
   if(mag.low>3 || mag.up<8) 
      data.select<-filter.mag(data.select,mag.low,mag.up)

   if(F.low>1.0 || F.up<3.0) 
      data.select<-filter.F(data.select,F.low,F.up)
 
   if(sol.low>0 || sol.up<360) 
      data.select<-filter.sol(data.select,sol.low,sol.up)

   if(h.low>10 || h.up<90) 
      data.select<-filter.h(data.select,shw,Ralpha,Delta,h.low,h.up)
   
   if(!is.null(r))
      data.select<-filter.totcor(data.select,shw,Ralpha,Delta,r,C)

   data.select
}
